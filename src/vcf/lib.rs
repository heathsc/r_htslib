use std::ptr::NonNull;
use std::convert::{AsRef, AsMut};

use super::{
   htsFile, kstring_t, HtsPos,
   BGZF, BcfRec, VcfHeader
};

use crate::bgzf_getline;

use c2rust_bitfields::BitfieldStruct;
use libc::{c_char, c_int, c_void};

pub const BCF_DT_ID: u32 = 0;
pub const BCF_DT_CTG: u32 = 1;
pub const BCF_DT_SAMPLE: u32 = 2;

#[repr(C)]
pub struct bcf_hdr_t {
   _unused: [u8; 0],
}

#[link(name = "hts")]
extern "C" {
   pub(super) fn bcf_hdr_init(mode: *const c_char) -> *mut bcf_hdr_t;
   pub(super) fn bcf_hdr_append(hdr: *mut bcf_hdr_t, line: *const c_char) -> c_int;
   pub(super) fn bcf_hdr_get_version(hdr: *const bcf_hdr_t) -> *const c_char;
   pub(super) fn bcf_hdr_add_sample(hdr: *mut bcf_hdr_t, sample: *const c_char) -> c_int;
   pub(super) fn bcf_hdr_write(fp: *mut htsFile, hdr: *mut bcf_hdr_t) -> c_int;
   pub(super) fn bcf_hdr_id2int(hdr: *const bcf_hdr_t, type_: c_int, id: *const c_char) -> c_int;
   pub(super) fn bcf_hdr_sync(hdr: *mut bcf_hdr_t) -> c_int;
   pub(super) fn bcf_init() -> *mut bcf1_t;
   pub(super) fn bcf_read(fp: *mut htsFile, hdr: *const bcf_hdr_t, v: *mut bcf1_t) -> c_int;
   pub(crate) fn bcf_readrec(fp: *mut BGZF, tbxv: *mut c_void, sv: *mut c_void, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos) -> c_int;
   pub(super) fn bcf_destroy(bcf: *mut bcf1_t);
   pub(super) fn bcf_clear(bcf: *mut bcf1_t);
   pub(super) fn bcf_write(hfile: *mut htsFile, hdr: *mut bcf_hdr_t, brec: *mut bcf1_t) -> c_int;
   pub(super) fn bcf_hdr_read(fp: *mut htsFile) -> *mut bcf_hdr_t;
   pub(super) fn bcf_hdr_destroy(hdr: *mut bcf_hdr_t);
   pub(super) fn vcf_format(hdr: *const bcf_hdr_t, v: *mut bcf1_t, s: *mut kstring_t) -> c_int;
   pub(super) fn vcf_parse(s: *mut kstring_t, hdr: *const bcf_hdr_t, v: *mut bcf1_t) -> c_int;
}

pub const BCF_BT_NULL: u8 = 0;
pub const BCF_BT_INT8: u8 = 1;
pub const BCF_BT_INT16: u8 = 2;
pub const BCF_BT_INT32: u8 = 3;
pub const BCF_BT_INT64: u8 = 4;
pub const BCF_BT_FLOAT: u8 = 5;
pub const BCF_BT_CHAR: u8 = 7;

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
struct bcf_fmt_t {
   id: c_int,
   n: c_int,
   size: c_int,
   _type: c_int,
   p: *mut u8,
   p_len: u32,
   #[bitfield(name = "p_off", ty = "u32", bits = "0..=30")]
   #[bitfield(name = "p_free", ty = "u8", bits = "31..=31")]
   bitfield1: [u8; 4],
}
#[repr(C)]
struct bcf_variant_t {
   _unused: [u8; 0],
}

#[repr(C)]
struct bcf_dec_t {
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
   fmt: *mut bcf_fmt_t,
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
   d: bcf_dec_t,
   max_unpack: c_int,
   unpacked: c_int,
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