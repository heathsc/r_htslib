use std::io;
use std::ptr::{NonNull, null, null_mut};
use std::ffi::{c_void, CStr, CString};
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

use lazy_static::lazy_static;
use regex::bytes::Regex;

use libc::{c_char, c_int};
use crate::bgzf_getline;
use super::{hts_err, hts_idx_t, HtsItr, Hts, HtsPos, HtsFileDesc, HtsRead, BGZF, kstring_t, hts_itr_next};

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

   pub fn parse(&mut self, tbx: &mut tbx_t) -> io::Result<()> {
      lazy_static! {
         static ref RE: Regex = Regex::new(r"^|;END=(\d+)").unwrap();
      }

      if let Some(p) = self.line.as_slice_mut(true) {
         let sc = tbx.conf.sc;
         let bc = tbx.conf.bc;
         let ec = tbx.conf.ec;
         let preset = tbx.conf.preset;
         let mut id = 0;
         let mut vals: (Option<c_int>, Option<HtsPos>, Option<HtsPos>) = (None, None, None);

         let get_ptr = |p: &[u8]| p.as_ptr() as *const u8 as *const c_char;

         let get_num = |p: &[u8]| -> io::Result<HtsPos> { unsafe {
            let p_start = get_ptr(p);
            let mut p_end: *mut c_char = null_mut();
            let k = libc::strtol(p_start, &mut p_end as *mut *mut c_char, 10) as HtsPos;
            if std::ptr::eq(p_start, p_end) {
               Err(hts_err("Error parsing integer".to_string()))
            } else {
               Ok(k)
            }
         }};

         // We can't just use split() here as we want to recover the delimiting character
         let mut p_rest = p;
         while !p_rest.is_empty() {
            let x = p_rest.iter().enumerate().find(|(_, &c)| c == b'\t').map(|(i, _)| i + 1).unwrap_or(p_rest.len());
            let (p_curr, q) = p_rest.split_at_mut(x);
            id += 1;
            let len =  p_curr.len();
            if len > 0 {
               if id == sc {
                  // Sequence (contig)
                  let mut tmp: u8 = 0;
                  vals.0 = unsafe {
                     std::ptr::swap(&mut p_curr[len - 1], &mut tmp);
                     match tbx_name2id(tbx as *mut tbx_t as *mut c_void, get_ptr(p_curr)) {
                        tid if tid >= 0 => Some(tid),
                        _ => return Err(hts_err("Out of memory".to_string())),
                     }
                  };
                  p_curr[len - 1] = tmp;

               } else if id == bc {
                  // Start column
                  let k = get_num(p_curr)?;
                  let (begin, end) = if (preset & TBX_UCSC) != 0 {
                     ((k - 1).max(0), k.max(1))
                  } else {
                     (k.max(0), (k + 1).max(1))
                  };
                  vals.1 = Some(begin);
                  vals.2 = Some(end);

               } else {
                  match preset & 0xffff {
                     TBX_GENERIC if id == ec => vals.2 = Some(get_num(p_curr)?),
                     TBX_SAM if id == 6 => { // CIGAR
                        let mut len = 0;
                        let mut op_len = 0;
                        for c in p_curr.iter() {
                           match *c {
                              x if x >= b'0' && x <= b'9' => op_len = (op_len * 10) + ((x - b'0') as HtsPos),
                              b'M' | b'D' | b'N' => {
                                 len += op_len;
                                 op_len = 0;
                              }
                              _ => op_len = 0,
                           }
                        }
                        vals.2 = vals.1.map(|k| k + len.max(1))
                     },
                     TBX_VCF if id == 4 => vals.2 = vals.1.map(|k| k + (len - 1) as HtsPos), // Use length of reference allele
                     TBX_VCF if id == 8 => { // Check for "END="
                        if let Some(r) = RE.captures(q).and_then(|cap| cap.get(1)).map(|x| x.as_bytes()) {
                           let end = get_num(r)?;
                           if let Some(begin) = vals.1 {
                              if begin < end { vals.2 = Some(end) }
                           }
                        }
                     },
                     _ => (),
                  }
               }
            }
            p_rest = q;
         }
         if let (Some(tid), Some(begin), Some(end)) = vals {
            self.tid = tid;
            self.begin = begin;
            self.end = end;
            Ok(())
         } else {
            Err(hts_err("Error parsing tbx record".to_string()))
         }
      } else {
         Err(hts_err("Can not parse tbx record from empty string".to_string()))
      }
   }

   pub fn to_str(&self) -> Option<&str> { self.line.to_str() }

   pub fn to_cstr(&self) -> Option<&CStr> { self.line.to_cstr() }

   pub fn tid(&self) -> usize { self.tid as usize }

   pub fn begin(&self) -> HtsPos { self.begin }

   pub fn end(&self) -> HtsPos { self.end }
}

impl HtsRead for TbxRec {
   fn read(&mut self, hts: &mut Hts) -> io::Result<bool> {
      let (fp, tbx) = hts.hts_file_and_tbx();
      if let Some(HtsFileDesc::Bgzf(bgzf)) = fp.file_desc() {
         let res = if let Some(tbx) = tbx {
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

   fn read_itr(&mut self, hts: &mut Hts, itr: &mut HtsItr) -> io::Result<bool> {
      let (fp, tbx) = hts.hts_file_and_tbx();
      if let Some(HtsFileDesc::Bgzf(bgzf)) = fp.file_desc() {
         let res = if let Some(tbx) = tbx {
            match unsafe { hts_itr_next(
               bgzf.as_ptr(),
               itr.as_mut(),
               self as *mut TbxRec as *mut c_void,
               tbx.as_mut() as *mut tbx_t as *mut c_void,
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
   pub(crate) fn tbx_index_load3(fname: *const c_char, fnidx: *const c_char, flags: c_int) -> *mut tbx_t;
   fn tbx_seqnames(tbx: *const tbx_t, n: *mut c_int) -> *mut *const c_char;
   pub (crate) fn tbx_readrec(fp: *mut BGZF, tbxv: *mut c_void, sv: *mut c_void, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos) -> c_int;
   fn tbx_destroy(tbx: *mut tbx_t);
   pub(crate) fn tbx_name2id(tbx: *mut c_void, ss: *const c_char) -> c_int;
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
   pub fn name2id(&mut self, s: &CStr) -> io::Result<usize> {
      let x = unsafe { tbx_name2id(self as *mut tbx_t as *mut c_void, s.as_ptr())};
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

pub (crate) unsafe extern "C" fn tbx_read_itr(fp: *mut BGZF, data: *mut c_void, r: *mut c_void, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos) -> c_int {
   let tbx = NonNull::new(data as *mut tbx_t).unwrap().as_mut();
   let tbx_rec = NonNull::new(r as *mut TbxRec).unwrap().as_mut();
   match bgzf_getline(fp, b'\n' as c_int, &mut tbx_rec.line) {
      c if c >= 0 => {
        match tbx_rec.parse(tbx) {
           Ok(_) => {
              *tid = tbx_rec.tid;
              *beg = tbx_rec.begin;
              *end = tbx_rec.end;
              c
           },
           Err(_) => -2,
        }
      },
      c => c,
   }
}
