use std::{
   ptr::NonNull,
   convert::{AsRef, AsMut, TryFrom},
   io
};

use super::{
   htsFile, kstring_t, HtsPos,
   BGZF, BcfRec, VcfHeader
};

use crate::{get_cstr, from_cstr, hts_err};

use crate::bgzf_getline;

use c2rust_bitfields::BitfieldStruct;
use libc::{c_char, c_int, c_void};

pub const BCF_DT_ID: u32 = 0;
pub const BCF_DT_CTG: u32 = 1;
pub const BCF_DT_SAMPLE: u32 = 2;

// Variable Length
pub const BCF_VL_FIXED: u32 = 0;
pub const BCF_VL_VAR: u32 = 1;
pub const BCF_VL_A: u32 = 2;
pub const BCF_VL_G: u32 = 3;
pub const BCF_VL_R: u32 = 4;

pub const BCF_UN_STR: usize = 1;    // up to ALT inclusive
pub const BCF_UN_FLT: usize = 2;    // up to FILTER
pub const BCF_UN_INFO: usize = 4;   // up to INFO
pub const BCF_UN_SHR: usize = BCF_UN_STR|BCF_UN_FLT|BCF_UN_INFO;    // All shared information
pub const BCF_UN_FMT: usize = 8;           // unpack format and each sample
pub const BCF_UN_IND: usize = BCF_UN_FMT;  // a synonym of BCF_UN_FMT
pub const BCF_UN_ALL: usize = BCF_UN_SHR|BCF_UN_FMT; // everything

pub const BCF_INT8_MISSING: i8 = i8::MIN;
pub const BCF_INT16_MISSING: i16 = i16::MIN;
pub const BCF_INT32_MISSING: i32 = i32::MIN;
pub const BCF_FLOAT_MISSING: u32 = 0x7f800001;

pub const BCF_INT8_VECTOR_END: i8 = i8::MIN + 1;
pub const BCF_INT16_VECTOR_END: i16 = i16::MIN + 1;
pub const BCF_INT32_VECTOR_END: i32 = i32::MIN + 1;
pub const BCF_FLOAT_VECTOR_END: u32 = 0x7f800002;

macro_rules! mk_try_from {
   ($en:ty, $base:ty, $t:ty) => {
      impl TryFrom<$t> for $en {
         type Error = &'static str;
         fn try_from(value: $t) -> Result<Self, Self::Error> {
            let value = <$base>::try_from(value).map_err(|_| "Error converting from $t to $en")?;
            Self::try_from(value)
         }
      }
   };

   ($en:ty, $base:ty, $t:ty, $($t1:ty),*) => {
      mk_try_from!($en, $base, $t);
      mk_try_from!($en, $base, $($t1),*);
   };
}

// Header Line

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum BcfHeaderLine {
   Flt = 0,
   Info = 1,
   Fmt = 2,
   Ctg = 3,
   Str = 4, // Structure header line TAG=<A..,B=..>
   Gen = 5, // Structure header line TAG=<A..,B=..>
}

impl TryFrom<u8> for BcfHeaderLine {
   type Error = &'static str;

   fn try_from(value: u8) -> Result<Self, Self::Error> {
      match value {
         0 => Ok(Self::Flt),
         1 => Ok(Self::Info),
         2 => Ok(Self::Fmt),
         3 => Ok(Self::Ctg),
         4 => Ok(Self::Str),
         5 => Ok(Self::Ctg),
         _ => Err("Bad BCF Header line"),
      }
   }
}

mk_try_from!(BcfHeaderLine, u8, i8, i16, i32, i64, isize, u16, u32, u64, usize);

// Header Line
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum BcfHeaderType {
   Flag = 0,
   Int = 1,
   Real= 2,
   Str = 3,
   Long = 0x101,
}

impl TryFrom<u16> for BcfHeaderType {
   type Error = &'static str;

   fn try_from(value: u16) -> Result<Self, Self::Error> {
      match value {
         0 => Ok(Self::Flag),
         1 => Ok(Self::Int),
         2 => Ok(Self::Real),
         3 => Ok(Self::Str),
         0x101 => Ok(Self::Long),
         _ => Err("Bad BCF Header type"),
      }
   }
}

mk_try_from!(BcfHeaderType, u16, i16, i32, i64, isize, u32, u64, usize);

pub enum BcfOpt<T> {
   Some(T),
   Missing,
   EndOfVec,
}

impl <T>BcfOpt<T> {
   pub fn is_some(&self) -> bool { matches!(self, Self::Some(_)) }

   pub fn is_missing(&self) -> bool { matches!(self, Self::Missing) }

   pub fn is_end_of_vec(&self) -> bool { matches!(self, Self::EndOfVec) }
}

fn try_parse_int(p: &[u8]) -> Result<BcfOpt<i32>, &'static str> {
   match p.len() {
      1 => {
         let x = i8::from_le_bytes([p[0]]);
         Ok(match x {
            BCF_INT8_MISSING => BcfOpt::Missing,
            BCF_INT8_VECTOR_END => BcfOpt::EndOfVec,
            _ => BcfOpt::Some(x as i32),
         })
      },
      2 => {
         let x = i16::from_le_bytes([p[0], p[1]]);
         Ok(match x {
            BCF_INT16_MISSING => BcfOpt::Missing,
            BCF_INT16_VECTOR_END => BcfOpt::EndOfVec,
            _ => BcfOpt::Some(x as i32),
         })
      },
      4 => {
         let x = i32::from_le_bytes([p[0], p[1], p[2], p[3]]);
         Ok(match x {
            BCF_INT32_MISSING => BcfOpt::Missing,
            BCF_INT32_VECTOR_END => BcfOpt::EndOfVec,
            _ => BcfOpt::Some(x),
         })
      },
      _ => Err("Illegal integer width"),
   }
}

pub trait BcfHeaderVar<'a> {
   type Item;

   fn hdr_type() -> BcfHeaderType;

   fn try_parse(p: &'a[u8]) -> Result<BcfOpt<Self::Item>, &'static str>;
}

pub struct BcfHdrFlag(bool);
impl <'a>BcfHeaderVar<'a> for BcfHdrFlag {
   type Item = bool;

   fn hdr_type() -> BcfHeaderType { BcfHeaderType::Flag }

   fn try_parse(p: &'a [u8]) -> Result<BcfOpt<Self::Item>, &'static str> {
      try_parse_int(p).map(|x| BcfOpt::Some(x.is_some()) )
   }
}

pub struct BcfHdrFloat(f32);
impl <'a>BcfHeaderVar<'a> for BcfHdrFloat {
   type Item = f32;

   fn hdr_type() -> BcfHeaderType { BcfHeaderType::Real }

   fn try_parse(p: &'a [u8]) -> Result<BcfOpt<Self::Item>, &'static str> {
      if p.len() == 4 {
         let x = u32::from_le_bytes([p[0], p[1], p[2], p[3]]);
         Ok(match x {
            BCF_FLOAT_MISSING => BcfOpt::Missing,
            BCF_FLOAT_VECTOR_END => BcfOpt::EndOfVec,
            _ => BcfOpt::Some(f32::from_le_bytes([p[0], p[1], p[2], p[3]])),
         })
      } else { Err("Illegal float width") }
   }
}

pub struct BcfHdrStr<'a>(&'a str);
impl <'a>BcfHeaderVar<'a> for BcfHdrStr<'a> {
   type Item = &'a str;

   fn hdr_type() -> BcfHeaderType { BcfHeaderType::Str }

   fn try_parse(p: &'a [u8]) -> Result<BcfOpt<Self::Item>, &'static str> {
      if p.is_empty() { Err("Empty string value") } else {
         let p = if let Some((ix, _)) = p.iter().enumerate().find(|(_, &c)| c == 0) {
            &p[..ix]
         } else { p };
         Ok(BcfOpt::Some(std::str::from_utf8(p).map_err(|_| "UTF8 coding not allowed in VCF/BCF strings")?))
      }
   }
}

pub struct BcfHdrInt(i32);
impl <'a>BcfHeaderVar<'a> for BcfHdrInt {
   type Item = i32;

   fn hdr_type() -> BcfHeaderType { BcfHeaderType::Int }

   fn try_parse(p: &'a [u8]) -> Result<BcfOpt<Self::Item>, &'static str> { try_parse_int(p) }
}

#[repr(C)]
pub struct bcf_hrec_t {
   _type: c_int,
   key: *mut c_char,
   val: *mut c_char,
   nkeys: c_int,
   keys: *mut *mut c_char,
   vals: *mut *mut c_char,
}

impl bcf_hrec_t {
   pub fn nkeys(&self) -> usize { self.nkeys as usize }
   pub fn find_key<S: AsRef<str>>(&self, key: S) -> Option<&str> {
      let ix = unsafe{bcf_hrec_find_key(self, get_cstr(key).as_ptr())};
      if ix < 0 {	None }
      else {
         assert!(ix < self.nkeys);
         Some(unsafe { from_cstr(*self.vals.add(ix as usize))})
      }
   }
   pub fn key(&self) -> &str { from_cstr(self.key) }
   pub fn val(&self) -> &str { from_cstr(self.val) }
   pub fn vals(&self, ix: usize) -> io::Result<&str> {
      if ix >= self.nkeys() { Err(hts_err("Invalid key id".to_string()))}
      else { Ok(from_cstr(unsafe{*self.vals.add(ix)})) }
   }
   pub fn keys(&self, ix: usize) -> io::Result<&str> {
      if ix >= self.nkeys() { Err(hts_err("Invalid key id".to_string()))}
      else { Ok(from_cstr(unsafe{*self.keys.add(ix)})) }
   }
   pub fn get_type(&self) -> c_int { self._type }
}

#[repr(C)]
pub struct bcf_idinfo_t {
   pub(super) info: [u64; 3],
   pub(super) hrec: [*mut bcf_hrec_t; 3],
   pub(super) id: c_int,
}

#[repr(C)]
pub struct bcf_idpair_t {
   pub(super) key: *const c_char,
   pub(super) val: *const bcf_idinfo_t,
}

#[repr(C)]
pub struct vdict_t {
   n_buckets: u32,
   size: u32,
   n_occupied: u32,
   upper_bound: u32,
   flags: *mut u32,
   keys: *mut c_void,
   vals: *mut c_void,
}

#[repr(C)]
pub struct bcf_hdr_t {
   pub(super) n: [i32; 3],
   pub(super) id: [*mut bcf_idpair_t; 3],
   pub(super) dict: [*mut vdict_t; 3],
   pub(super) samples: *mut *mut c_char,
   pub(super) hrec: *mut *mut bcf_hrec_t,
   pub(super) nhrec: c_int,
   dirty: c_int,
   ntransl: c_int,
   transl: [*mut c_int; 2],
   n_samples_ori: c_int,
   keep_samples: *mut u8,
   mem: kstring_t,
   m: [i32; 3],
}

impl bcf_hdr_t {
   pub fn n_samples(&self) -> i32 { self.n[BCF_DT_SAMPLE as usize] }
}

#[link(name = "hts")]
extern "C" {
   pub(super) fn bcf_hdr_init(mode: *const c_char) -> *mut bcf_hdr_t;
   pub(super) fn bcf_hdr_append(hdr: *mut bcf_hdr_t, line: *const c_char) -> c_int;
   pub(super) fn bcf_hdr_dup(hdr: *const bcf_hdr_t) -> *mut bcf_hdr_t;
   pub(super) fn bcf_hdr_get_version(hdr: *const bcf_hdr_t) -> *const c_char;
   pub(super) fn bcf_hdr_add_sample(hdr: *mut bcf_hdr_t, sample: *const c_char) -> c_int;
   pub(super) fn bcf_hdr_write(fp: *mut htsFile, hdr: *mut bcf_hdr_t) -> c_int;
   pub(super) fn bcf_hdr_id2int(hdr: *const bcf_hdr_t, type_: c_int, id: *const c_char) -> c_int;
   pub(super) fn bcf_hdr_sync(hdr: *mut bcf_hdr_t) -> c_int;
   pub(super) fn bcf_hdr_seqnames(hdr: *const bcf_hdr_t, n_seqs: *mut c_int) -> *mut *const c_char;
   pub(super) fn bcf_init() -> *mut bcf1_t;
   pub(super) fn bcf_read(fp: *mut htsFile, hdr: *const bcf_hdr_t, v: *mut bcf1_t) -> c_int;
   pub(crate) fn bcf_readrec(fp: *mut BGZF, tbxv: *mut c_void, sv: *mut c_void, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos) -> c_int;
   pub(super) fn bcf_destroy(bcf: *mut bcf1_t);
   pub(super) fn bcf_clear(bcf: *mut bcf1_t);
   pub(super) fn bcf_unpack(b: *mut bcf1_t, which: c_int);
   pub(super) fn bcf_write(hfile: *mut htsFile, hdr: *mut bcf_hdr_t, brec: *mut bcf1_t) -> c_int;
   pub(super) fn bcf_hdr_read(fp: *mut htsFile) -> *mut bcf_hdr_t;
   pub(super) fn bcf_hdr_destroy(hdr: *mut bcf_hdr_t);
   pub(super) fn bcf_hrec_find_key(hrec: *const bcf_hrec_t, key: *const c_char) -> c_int;
   pub(super) fn vcf_format(hdr: *const bcf_hdr_t, v: *mut bcf1_t, s: *mut kstring_t) -> c_int;
   pub(super) fn vcf_parse(s: *mut kstring_t, hdr: *const bcf_hdr_t, v: *mut bcf1_t) -> c_int;
}

pub const BCF_BT_NULL: c_int = 0;
pub const BCF_BT_INT8: c_int = 1;
pub const BCF_BT_INT16: c_int = 2;
pub const BCF_BT_INT32: c_int = 3;
pub const BCF_BT_INT64: c_int = 4;
pub const BCF_BT_FLOAT: c_int = 5;
pub const BCF_BT_CHAR: c_int = 7;

pub const BCF_MAX_BT_INT8: i32 = 0x7f; /* INT8_MAX  */
pub const BCF_MAX_BT_INT16: i32 = 0x7fff; /* INT16_MAX */
pub const MAX_BT_INT32: i32 = 0x7fffffff; /* INT32_MAX */
pub const BCF_MIN_BT_INT8: i32 = -120; /* INT8_MIN  + 8 */
pub const BCF_MIN_BT_INT16: i32 = -32760; /* INT16_MIN + 8 */
pub const BCF_MIN_BT_INT32: i32 = -2147483640; /* INT32_MIN + 8 */

#[allow(non_upper_case_globals)]
pub const bcf_int8_vector_end: i32 = -127; /* INT8_MIN  + 1 */
#[allow(non_upper_case_globals)]
pub const bcf_int16_vector_end: i32 = -32767; /* INT16_MIN + 1 */
#[allow(non_upper_case_globals)]
pub const bcf_int32_vector_end: i32 = -2147483647; /* INT32_MIN + 1 */
#[allow(non_upper_case_globals)]
pub const bcf_int64_vector_end: i64 = -9223372036854775807; /* INT64_MIN + 1 */
#[allow(non_upper_case_globals)]
pub const bcf_str_vector_end: usize = 0;
#[allow(non_upper_case_globals)]
pub const bcf_int8_missing: i32 = -128; /* INT8_MIN  */
#[allow(non_upper_case_globals)]
pub const bcf_int16_missing: i32 = -32767 - 1; /* INT16_MIN */
#[allow(non_upper_case_globals)]
pub const bcf_int32_missing: i32 = -2147483647 - 1; /* INT32_MIN */
#[allow(non_upper_case_globals)]
pub const bcf_int64_missing: i64 = -9223372036854775807 - 1; /* INT64_MIN */
#[allow(non_upper_case_globals)]
pub const bcf_str_missing: usize = 0x07;

#[repr(C)]
struct bcf_info_t {
   _unused: [u8; 0],
}
#[repr(C)]
#[derive(BitfieldStruct)]
pub(super) struct bcf_fmt_t {
   pub(super) id: c_int,
   pub(super) n: c_int,
   pub(super) size: c_int,
   pub(super) vtype: c_int,
   pub(super) p: *mut u8,
   pub(super) p_len: u32,
   #[bitfield(name = "p_off", ty = "u32", bits = "0..=30")]
   #[bitfield(name = "p_free", ty = "u8", bits = "31..=31")]
   pub(super) bitfield1: [u8; 4],
}

impl bcf_fmt_t {
   pub(super) fn type_size(&self) -> usize {
      match self.vtype {
         BCF_BT_INT8 | BCF_BT_CHAR => 1,
         BCF_BT_INT16 => 2,
         BCF_BT_INT32 | BCF_BT_FLOAT => 4,
         BCF_BT_INT64 => 8,
         _ => 0,
      }
   }
}

#[repr(C)]
struct bcf_variant_t {
   _unused: [u8; 0],
}

#[repr(C)]
pub(super) struct bcf_dec_t {
   m_fmt: c_int,
   m_info: c_int,
   m_id: c_int,
   m_als: c_int,
   m_allele: c_int,
   m_flt: c_int,
   n_flt: c_int,
   flt: *mut c_int,
   id: *mut c_char,
   als: *mut c_char,
   alleles: *mut *mut c_char,
   info: *mut bcf_info_t,
   pub(super) fmt: *mut bcf_fmt_t,
   var: *mut bcf_variant_t,
   n_var: c_int,
   var_type: c_int,
   shared_dirty: c_int,
   indiv_dirty: c_int,
}

#[repr(C)]
#[derive(BitfieldStruct)]
pub struct bcf1_t {
   pos: HtsPos,
   rlen: HtsPos,
   rid: i32,
   pub(super) qual: f32,
   #[bitfield(name = "n_info", ty = "u16", bits = "0..=15")]
   #[bitfield(name = "n_allele", ty = "u16", bits = "16..=31")]
   bitfield1: [u8; 4],
   #[bitfield(name = "n_fmt", ty = "u8", bits = "0..=7")]
   #[bitfield(name = "n_sample", ty = "u32", bits = "8..=31")]
   bitfield2: [u8; 4],
   shared: kstring_t,
   indiv: kstring_t,
   pub(super) d: bcf_dec_t,
   max_unpack: c_int,
   pub(super) unpacked: c_int,
   unpack_size: [c_int; 3],
   pub(super) errcode: c_int,
}

impl bcf1_t {
   pub fn clear(&mut self) {
      unsafe { bcf_clear(self) }
   }

   pub fn shared(&mut self) -> &mut kstring_t {
      &mut self.shared
   }
   pub fn indiv(&mut self) -> &mut kstring_t {
      &mut self.indiv
   }
   pub fn set_rid(&mut self, rid: usize) {
      self.rid = rid as i32
   }
   pub fn set_pos(&mut self, pos: usize) {
      self.pos = pos as HtsPos
   }
   pub fn set_rlen(&mut self, rlen: usize) {
      self.rlen = rlen as HtsPos
   }
   pub fn set_qual(&mut self, qual: f32) {
      self.qual = qual
   }
   pub fn format(&mut self, hdr: &VcfHeader, s: &mut kstring_t) -> c_int {
      unsafe { vcf_format(hdr.as_ref(), self, s) }
   }

   pub fn unpack(&mut self, which: usize) { unsafe{bcf_unpack(self, which as c_int)} }
}

pub (crate) unsafe extern "C" fn vcf_read_itr(fp: *mut BGZF, data: *mut c_void, r: *mut c_void, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos) -> c_int {
   let hdr = NonNull::new(data as *mut bcf_hdr_t).unwrap().as_ref();
   let bcf_rec = NonNull::new(r as *mut BcfRec).unwrap().as_mut();
   match bgzf_getline(fp, b'\n' as c_int, &mut bcf_rec.line) {
      c if c >= 0 => {
         match vcf_parse(&mut bcf_rec.line, hdr, bcf_rec.as_mut()) {
            0 => {
               *tid = bcf_rec.rid;
               *beg = bcf_rec.pos;
               *end = bcf_rec.pos + bcf_rec.rlen;
               c
            },
            c => c,
         }
      },
      c => c,
   }
}

pub (crate) unsafe extern "C" fn bcf_read_itr(fp: *mut BGZF, data: *mut c_void, r: *mut c_void, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos) -> c_int {
   let bcf_rec = NonNull::new(r as *mut BcfRec).unwrap().as_mut();
   bcf_readrec(fp, data, bcf_rec.as_mut() as *mut bcf1_t as *mut c_void, tid, beg, end)
}

pub (crate) unsafe extern "C" fn bcf_hdr_name2id(hdr: *mut c_void, s: *const c_char) -> c_int {
   bcf_hdr_id2int(hdr as *const bcf_hdr_t, BCF_DT_CTG as c_int, s)
}