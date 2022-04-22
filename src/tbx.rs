use std::io;
use std::ptr::{NonNull, null};
use std::ffi::{c_void, CStr, CString};
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

use libc::{c_char, c_int};
use crate::HtsItrType;
use super::{hts_err, hts_idx_t, hts_itr_t, HtsItr, HtsPos, HtsHdr, HtsFile, HtsFileDesc, HtsRead, BGZF, kstring_t};

pub const TBX_GENERIC: i32 = 0;
pub const TBX_SAM: i32 = 1;
pub const TBX_VCF: i32 = 2;
pub const TBX_UCSC: i32 = 0x10000;

#[repr(C)]
pub struct tbx_conf_t {
   preset: i32,
   sc: i32,
   bc: i32,
   ec: i32,
   meta_char: i32,
   line_skip: i32,
}

#[allow(non_upper_case_globals)]
pub const tbx_conf_vcf: tbx_conf_t = tbx_conf_t{ preset: TBX_VCF, sc: 1, bc: 2, ec: 0, meta_char: b'#' as i32, line_skip: 0};

#[repr(C)]
pub struct tbx_t {
   conf: tbx_conf_t,
   index: NonNull<hts_idx_t>,
   dict: *mut c_void,
}

#[derive(Debug, Default)]
#[repr(C)]
pub struct TbxRec {
   tid: c_int,
   begin: HtsPos,
   end: HtsPos,
   line: kstring_t,
}

unsafe impl Send for TbxRec {}

impl TbxRec {
   pub fn new() -> Self { Self::default() }

   pub fn to_str(&self) -> Option<&str> { self.line.to_str() }

   pub fn to_cstr(&self) -> Option<&CStr> { self.line.to_cstr() }

   pub fn tid(&self) -> usize { self.tid as usize }

   pub fn begin(&self) -> HtsPos { self.begin }

   pub fn end(&self) -> HtsPos { self.end }
}

impl HtsRead for TbxRec {
   fn read(&mut self, fp: &mut HtsFile, hdr: Option<&mut HtsHdr>) -> io::Result<bool> {
      if let Some(HtsFileDesc::Bgzf(bgzf)) = fp.file_desc() {
         let res = if let Some(HtsHdr::Tbx(tbx)) = hdr {
            match unsafe { tbx_readrec(
               bgzf.as_ptr(),
               tbx.as_mut() as *mut tbx_t as *mut c_void,
               &mut self.line as *mut kstring_t as *mut c_void,
               &mut self.tid as *mut c_int,
               &mut self.begin as *mut HtsPos,
               &mut self.end as *mut HtsPos
            )} {
               0..=c_int::MAX => Ok(true),
               -1 => Ok(false),
               _ => Err(hts_err("Error reading SAM/BAM record".to_string())),
            }
         } else {
            Err(hts_err(format!("Appropriate header not found for file {}", fp.name())))
         };
         res
      } else {
         Err(hts_err(format!("File {} is not in bgzf format (required for tabix)", fp.name())))
      }
   }
}

#[link(name = "hts")]
extern "C" {
   fn hts_itr_query(idx: *const hts_idx_t, tid: c_int, beg: HtsPos, end: HtsPos,
                    readrec: unsafe extern "C" fn (*mut BGZF, *mut c_void, *mut c_void, *mut c_int, *mut HtsPos, *mut HtsPos) -> c_int) -> *mut hts_itr_t;
   pub(crate) fn tbx_index_load3(fname: *const c_char, fnidx: *const c_char, flags: c_int) -> *mut tbx_t;
   fn tbx_seqnames(tbx: *const tbx_t, n: *mut c_int) -> *mut *const c_char;
   fn tbx_readrec(fp: *mut BGZF, tbxv: *mut c_void, sv: *mut c_void, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos) -> c_int;
   fn tbx_destroy(tbx: *mut tbx_t);
   fn tbx_name2id(tbx: *const tbx_t, ss: *const c_char) -> c_int;
   pub fn tbx_index_build(fname: *const c_char, min_shift: c_int, conf: *const tbx_conf_t) -> c_int;
}

impl tbx_t {
   pub fn seq_names(&self) -> Option<Vec<&str>> {
      let mut n_seq: c_int = 0;
      let p = unsafe{tbx_seqnames(self, &mut n_seq as *mut c_int)};
      if p.is_null() {
         None
      } else {
         let mut v = Vec::with_capacity(n_seq as usize);
         for i in 0..n_seq {
            let c_str: &CStr = unsafe { CStr::from_ptr(*p.offset(i as isize)) };
            let str_slice: &str = c_str.to_str().unwrap();
            v.push(str_slice);
         }
         unsafe {libc::free(p as *mut c_void)};
         Some(v)
      }
   }
   pub fn tbx_itr_queryi(&mut self, tid: c_int, beg: HtsPos, end: HtsPos) -> io::Result<HtsItr> {
      HtsItr::new(unsafe {hts_itr_query(self.index.as_ref(), tid, beg, end, tbx_readrec)}, HtsItrType::TbxItr(self)).ok_or_else(|| hts_err("Failed to obtain tbx iterator".to_string()))
   }
   pub fn name2id(&self, s: &CStr) -> io::Result<usize> {
      let x = unsafe { tbx_name2id(self, s.as_ptr())};
      if x >= 0 {
         Ok(x as usize)
      } else {
         Err(hts_err(format!("Failed to convert {:?} to id", s)))
      }
   }
   pub fn idx(&self) -> &hts_idx_t { unsafe { self.index.as_ref() } }
}

#[derive(Debug)]
pub struct Tbx {
   inner: NonNull<tbx_t>,
   phantom: PhantomData<tbx_t>,
}

impl Deref for Tbx {
   type Target = tbx_t;
   #[inline]
   fn deref(&self) -> &tbx_t { unsafe{self.inner.as_ref()} }
}

impl DerefMut for Tbx {
   #[inline]
   fn deref_mut(&mut self) -> &mut tbx_t {unsafe{ self.inner.as_mut() }}
}

impl AsRef<tbx_t> for Tbx {
   #[inline]
   fn as_ref(&self) -> &tbx_t { self}
}

impl AsMut<tbx_t> for Tbx {
   #[inline]
   fn as_mut(&mut self) -> &mut tbx_t { self}
}

impl Drop for Tbx {
   fn drop(&mut self) {
      unsafe { tbx_destroy(self.as_mut()) };
   }
}

impl Tbx {
   pub fn new(inner: NonNull<tbx_t>) -> Self { Self { inner, phantom: PhantomData }}

   pub fn load<S: AsRef<Path>>(name: S) -> io::Result<Self> {
      let name = name.as_ref();
      let cname = CString::new(name.as_os_str().as_bytes()).unwrap();
      match NonNull::new(unsafe{ tbx_index_load3(cname.as_ptr(), null::<c_char>(), 0)}) {
         None =>	Err(hts_err(format!("Couldn't open tabix index for file {}", name.display()))),
         Some(p) => Ok(Tbx{inner: p, phantom: PhantomData}),
      }
   }
}

pub enum TbxReadResult {
   Ok,
   EOF,
   Error,
}
