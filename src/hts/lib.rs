use c2rust_bitfields::BitfieldStruct;
use libc::{c_char, c_int, c_short, c_uchar, c_uint, off_t, size_t, ssize_t};
use std::ffi::{c_void};
use std::io;
use std::marker::PhantomData;
use std::ptr::{NonNull};


use crate::{
   HtsThreadPool,
   get_cstr, from_cstr, hts_err, kstring_t, sam_hdr_t,
};

pub type HtsPos = i64;

pub const HTS_IDX_NOCOOR: c_int = -2;
pub const HTS_IDX_START: c_int = -3;
pub const HTS_IDX_REST: c_int = -4;
pub const HTS_IDX_NONE: c_int = -5;

pub const HTS_IDX_SILENT_FAIL: c_int = 2;

pub const HTS_FMT_CSI: c_int = 0;
pub const HTS_FMT_BAI: c_int = 1;
pub const HTS_FMT_TBI: c_int = 2;
pub const HTS_FMT_CRAI: c_int = 3;
pub const HTS_FMT_FAI: c_int = 4;

pub const FT_UNKN: u32 = 0;
pub const FT_GZ: u32 = 1;
pub const FT_VCF: u32 = 2;
pub const FT_VCF_GZ: u32 = FT_GZ | FT_VCF;
pub const FT_BCF: u32 = 4;
pub const FT_BCF_GZ: u32 = FT_GZ | FT_BCF;
pub const FT_STDIN: u32 = 8;

// Region Parsing
pub const HTS_PARSE_THOUSANDS_SEP: c_int = 1; // Ignore ',' separators within numbers
pub const HTS_PARSE_ONE_COORD: c_int = 2; // chr:pos = chr:pos-pos not chr:pos-end
pub const HTS_PARSE_LIST: c_int = 4; // Expects a comma separated list of regions (disables HTS_PARSE_THOUSANDS_SEP)

#[repr(C)]
pub struct hts_opt {
   _unused: [u8; 0],
}

#[repr(C)]
pub struct hts_filter_t {
   _unused: [u8; 0],
}

#[repr(C)]
pub(super) struct hts_tpool {
   _unused: [u8; 0],
}

#[repr(C)]
pub struct prehtsFile {
   bfield: [u8; 4],
   lineno: i64,
   line: kstring_t,
   pub(super) fn_: *mut c_char,
}

#[repr(C)]
pub union FileType {
   pub(crate) bgzf: *mut BGZF,
   pub(crate) cram_fd: *mut cram_fd,
   pub(crate) hfile: *mut hfile,
}

#[repr(C)]
#[derive(BitfieldStruct)]
pub struct htsFile {
   #[bitfield(name = "is_bin", ty = "c_uchar", bits = "0..=0")]
   #[bitfield(name = "is_write", ty = "c_uchar", bits = "1..=1")]
   #[bitfield(name = "is_be", ty = "c_uchar", bits = "2..=2")]
   #[bitfield(name = "is_cram", ty = "c_uchar", bits = "3..=3")]
   #[bitfield(name = "is_bgzf", ty = "c_uchar", bits = "4..=4")]
   #[bitfield(name = "dummy", ty = "u32", bits = "5..=31")]
   bfield: [u8; 4],
   lineno: i64,
   line: kstring_t,
   fn_: NonNull<c_char>,
   fn_aux: *mut c_char,
   pub (crate) fp: FileType,
   state: *mut c_void,
   format: htsFormat,
   idx: *mut hts_idx_t,
   fnidx: *const c_char,
   bam_header: *mut sam_hdr_t,
   filter: *mut hts_filter_t,
}

pub enum HtsFileDesc {
   Bgzf(NonNull<BGZF>),
   Cram(NonNull<cram_fd>),
   Hfile(NonNull<hfile>),
}

impl htsFile {
   pub fn format(&self) -> &htsFormat { &self.format }

   pub fn name(&self) -> &str {
      from_cstr(self.name_ptr())
   }

   pub fn file_desc(&mut self) -> Option<HtsFileDesc> {
      unsafe { if self.is_bgzf() != 0 {
         NonNull::new(self.fp.bgzf).map(|p| HtsFileDesc::Bgzf(p))
      } else if self.is_cram() != 0 {
         NonNull::new(self.fp.cram_fd).map(|p| HtsFileDesc::Cram(p))
      } else {
         NonNull::new(self.fp.hfile).map(|p| HtsFileDesc::Hfile(p))
      }}
   }

   pub (super) fn name_ptr(&self) -> *const c_char {
      self.fn_.as_ptr()
   }

   pub fn set_threads(&mut self, t: usize) -> io::Result<()> {
      let ret = unsafe { hts_set_threads(self, t as c_int) };
      if ret != 0 {
         Err(hts_err(format!("Failed to set threads for file {}", self.name())))
      } else {
         Ok(())
      }
   }

   pub fn set_fai_filename<S: AsRef<str>>(&mut self, name: S) -> io::Result<()> {
      let name = name.as_ref();
      let ret = unsafe { hts_set_fai_filename(self, get_cstr(name).as_ptr()) };
      if ret != 0 {
         Err(hts_err(format!("Failed to attach reference index {} to file {}", name, self.name())))
      } else {
         Ok(())
      }
   }
}

#[repr(C)]
pub struct hts_idx_t {
   _unused: [u8; 0],
}

#[repr(C)]
pub struct hts_reglist_t {
   _unused: [u8; 0],
}

#[repr(C)]
pub struct hts_pair64_max {
   _unused: [u8; 0],
}

#[repr(C)]
pub struct HtsItrBins {
   n: c_int,
   m: c_int,
   a: *mut c_int,
}

#[repr(C)]
#[derive(BitfieldStruct)]
pub struct hts_itr_t {
   #[bitfield(name = "read_rest", ty = "c_uchar", bits = "0..=0")]
   #[bitfield(name = "finished", ty = "c_uchar", bits = "1..=1")]
   #[bitfield(name = "is_cram", ty = "c_uchar", bits = "2..=2")]
   #[bitfield(name = "nocoor", ty = "c_uchar", bits = "3..=3")]
   #[bitfield(name = "multi", ty = "c_uchar", bits = "4..=4")]
   #[bitfield(name = "dummy", ty = "u32", bits = "5..=31")]
   bfield: [u8; 4],
   tid: c_int,
   n_off: c_int,
   i: c_int,
   n_reg: c_int,
   beg: HtsPos,
   end: HtsPos,
   reg_list: *mut hts_reglist_t,
   curr_tid: c_int,
   curr_reg: c_int,
   curr_intv: c_int,
   curr_beg: HtsPos,
   curr_end: HtsPos,
   curr_off: u64,
   nocoor_off: u64,
   off: *mut hts_pair64_max,
   readrec: HtsReadrecFunc,
   seek: HtsSeekFunc,
   tell: HtsTellFunc,
   bins: HtsItrBins,
}

pub type HtsReadrecFunc = unsafe extern fn(fp: *mut BGZF, data: *mut c_void, r: *mut c_void, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos) -> c_int;
pub type HtsSeekFunc = unsafe extern fn(fp: *mut c_void, offset: i64, where_: c_int);
pub type HtsTellFunc = unsafe extern fn(fp: *mut c_void);
pub type HtsItrQueryFunc = unsafe extern fn(idx: *const hts_idx_t, tid: c_int, beg: HtsPos, end: HtsPos, readrec: HtsReadrecFunc) -> c_int;
pub type HtsName2Id = unsafe extern fn(hdr: *mut c_void, str: *const c_char) -> c_int;

#[repr(C)]
pub struct BGZF {
   _unused: [u8; 0],
}

#[repr(C)]
pub struct cram_fd {
   _unused: [u8; 0],
}

#[repr(C)]
pub struct hFile_backend {
   _unused: [u8; 0],
}

const UINT_LEN: usize = std::mem::size_of::<c_uint>();

#[repr(C)]
#[derive(BitfieldStruct)]
pub struct hfile {
   buffer: *mut c_char,
   begin: *mut c_char,
   end: *mut c_char,
   limit: *mut c_char,
   backend: *mut hFile_backend,
   offset: off_t,
   #[bitfield(name = "at_eof", ty = "c_uint", bits = "0..=0")]
   #[bitfield(name = "mobile", ty = "c_uint", bits = "1..=1")]
   #[bitfield(name = "readonly", ty = "c_uint", bits = "2..=2")]
   bfield: [u8; UINT_LEN],
   has_errno: c_int,
}

#[repr(C)]
pub enum htsFormatCategory {
   UnknownCategory,
   SequenceData,
   VariantData,
   IndexFile,
   RegionList,
}

impl Default for htsFormatCategory {
   fn default() -> Self {
      Self::UnknownCategory
   }
}

#[repr(C)]
#[derive(PartialEq)]
pub enum htsExactFormat {
   UnknownFormat, BinaryFormat, TextFormat,
   Sam, Bam, Bai, Cram, Crai, Vcf, Bcf, Csi, Gzi, Tbi, Bed, HtsGet,
   EmptyFormat, FastaFormat, FastqFormat, FaiFormat, FqiFormat, HtsCrypt4GH,
}

impl Default for htsExactFormat {
   fn default() -> Self { Self::UnknownFormat }
}

#[repr(C)]
#[derive(PartialEq)]
pub enum htsCompression {
   NoCompression,
   Gzip,
   Bgzf,
   Custom,
   Bzip2Compression,
}

impl Default for htsCompression {
   fn default() -> Self {
      Self::NoCompression
   }
}

#[repr(C)]
#[derive(Default)]
pub struct htsFormatVersion {
   major: c_short,
   minor: c_short,
}

#[repr(C)]
pub struct htsFormat {
   category: htsFormatCategory,
   pub(super) format: htsExactFormat,
   version: htsFormatVersion,
   pub(super) compression: htsCompression,
   compression_level: c_short,
   pub(super) specific: *mut hts_opt,
   phantom: PhantomData<hts_opt>,
}

impl Default for htsFormat {
   fn default() -> Self {
      Self {
         category: htsFormatCategory::default(),
         format: htsExactFormat::default(),
         version: htsFormatVersion::default(),
         compression: htsCompression::default(),
         compression_level: 0,
         specific: std::ptr::null_mut::<hts_opt>(),
         phantom: PhantomData,
      }
   }
}

impl Drop for htsFormat {
   fn drop(&mut self) {
      if !self.specific.is_null() {
         unsafe { hts_opt_free(self.specific) }
      }
   }
}

impl htsFormat {
   pub fn format(&self) -> &htsExactFormat { &self.format }
}

#[link(name = "hts")]
extern "C" {
   pub(super) fn hts_open_format(fn_: *const c_char, mode: *const c_char, fmt: *const htsFormat) -> *mut prehtsFile;
   pub(super) fn hts_close(fp_: *mut htsFile) -> c_int;
   pub(super) fn hts_set_threads(fp_: *mut htsFile, t_: c_int) -> c_int;
   pub(super) fn hts_set_thread_pool(fp_: *mut htsFile, p: *const HtsThreadPool) -> c_int;
   pub(super) fn hts_tpool_init(n: c_int) -> *mut hts_tpool;
   pub(super) fn hts_tpool_destroy(p: *mut hts_tpool);
   pub(super) fn sam_index_load3(fp_: *mut htsFile, name: *const c_char, fnidx: *const c_char, flags: c_int) -> *mut hts_idx_t;
   pub(super) fn hts_idx_load3(name: *const c_char, fnidx: *const c_char, fmt: c_int, flags : c_int) -> *mut hts_idx_t;
   pub(super) fn hts_idx_destroy(idx: *mut hts_idx_t);
   pub(super) fn hts_set_fai_filename(fp_: *mut htsFile, fn_aux: *const c_char) -> c_int;
   pub(super) fn hts_itr_destroy(iter: *mut hts_itr_t);
   pub(super) fn hts_itr_query(idx: *const hts_idx_t, tid: c_int, begin: HtsPos, end: HtsPos, readrec: HtsReadrecFunc) -> *mut hts_itr_t;
   pub(super) fn hts_parse_region(s: *const c_char, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos, getid: HtsName2Id, hdr: *mut c_void, flags: c_int) -> *const c_char;
   pub(super) fn hts_parse_format(format: *mut htsFormat, str: *const c_char) -> c_int;
   pub(super) fn hts_opt_free(opts: *mut hts_opt);
   pub(super) fn hts_opt_add(opts: *mut *mut hts_opt, c_arg: *const c_char) -> c_int;
   pub(crate) fn bgzf_getline(fp: *mut BGZF, delim: c_int, str: *mut kstring_t) -> c_int;
   pub(crate) fn bgzf_write(fp: *mut BGZF, data: *const c_void, len: size_t) -> ssize_t;
   pub(crate) fn bgzf_mt(fp: *mut BGZF, n_threads: c_int, n_sub_blocks: c_int) -> c_int;
//   pub(crate) fn bgzf_block_write(fp: *mut BGZF, data: *const c_void, len: size_t) -> ssize_t;
   pub(crate) fn bgzf_flush(fp: *mut BGZF) -> c_int;
   pub(crate) fn hflush(fp: &mut hfile) -> c_int;
   fn hfile_set_blksize(fp: *mut hfile, bufsize: size_t) -> c_int;
   fn hwrite2(fp: *mut hfile, data: *const c_void, total: size_t, copied: size_t) -> ssize_t;
   pub(crate) fn hts_itr_multi_next(fp: *mut htsFile, itr: *mut hts_itr_t, r: *mut c_void) -> c_int;
   pub(crate) fn hts_itr_next(fp: *mut BGZF, itr: *mut hts_itr_t, r: *mut c_void, data: *mut c_void) -> c_int;
}

pub fn hwrite(fp: &mut hfile, data: *const u8, nbytes: size_t) -> ssize_t {
   let nbytes1 = nbytes as isize;
   if fp.mobile() == 0 {
      let n = unsafe { fp.limit.offset_from(fp.begin) };
      if n < nbytes1 {
         unsafe { hfile_set_blksize(fp, (fp.limit.offset_from(fp.buffer) + nbytes1) as size_t ); }
         fp.end = fp.limit;
      }
   }
   let n = unsafe { fp.limit.offset_from(fp.begin) };
   if nbytes1 >= n && fp.begin == fp.buffer {
      // Go straight to hwrite2 if the buffer is empty and the request won't fit
      unsafe { hwrite2(fp, data as *const c_void, nbytes, 0) }
   } else {
      let n = n.min(nbytes1) as size_t;
      unsafe {
         libc::memcpy(fp.begin as *mut c_void, data as *const c_void, n);
         fp.begin = fp.begin.offset(n as isize)
      }
      if n == nbytes {
         n as ssize_t
      } else {
         unsafe {hwrite2(fp, data as *const c_void, nbytes, n) }
      }
   }

}
