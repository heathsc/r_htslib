use std::convert::TryInto;
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};
use std::ptr::{NonNull, null_mut};

use std::{fmt, io};

use crate::{
   get_cstr, htsFile, hts_err, kstring_t, Hts, HtsFile,
   HtsHdr, HtsItr, HtsRead, BGZF, hts_itr_next, hts_itr_multi_next,
};

use libc::{c_char, c_int, c_void};
pub mod lib;
pub use lib::*;

#[derive(Debug)]
pub struct SamHeader {
   inner: NonNull<sam_hdr_t>,
   phantom: PhantomData<sam_hdr_t>,
}

impl Deref for SamHeader {
   type Target = sam_hdr_t;
   #[inline]
   fn deref(&self) -> &sam_hdr_t { unsafe{self.inner.as_ref()} }
}

impl DerefMut for SamHeader {
   #[inline]
   fn deref_mut(&mut self) -> &mut sam_hdr_t {unsafe{ self.inner.as_mut() }}
}

impl AsRef<sam_hdr_t> for SamHeader {
   #[inline]
   fn as_ref(&self) -> &sam_hdr_t { self}
}

impl AsMut<sam_hdr_t> for SamHeader {
   #[inline]
   fn as_mut(&mut self) -> &mut sam_hdr_t { self}
}

impl Drop for SamHeader {
   fn drop(&mut self) {
      unsafe { sam_hdr_destroy(self.as_mut()) };
   }
}

impl Default for SamHeader {
   fn default() -> Self {
      let inner =
         NonNull::new(unsafe { sam_hdr_init() }).expect("Failed to allocated new SAM header");
      Self {
         inner,
         phantom: PhantomData,
      }
   }
}

impl SamHeader {
   pub fn new() -> Self {
      Self::default()
   }

   pub fn read(hts_file: &mut HtsFile) -> io::Result<Self> {
      match NonNull::new(unsafe { sam_hdr_read(hts_file.as_mut()) }) {
         None => Err(hts_err(format!("Failed to load header from {}", hts_file.name()))),
         Some(p) => Ok(Self {
            inner: p,
            phantom: PhantomData,
         }),
      }
   }

   pub fn dup(&self) -> io::Result<Self> {
      match NonNull::new(unsafe { sam_hdr_dup(self.as_ref()) }) {
         None => Err(hts_err("Failed to duplicate SAM/BAM header".to_string())),
         Some(p) => Ok(SamHeader {
            inner: p,
            phantom: PhantomData,
         }),
      }
   }
}


#[derive(Debug)]
pub struct BamRec {
   inner: NonNull<bam1_t>,
   phantom: PhantomData<bam1_t>,
}

impl Deref for BamRec {
   type Target = bam1_t;
   #[inline]
   fn deref(&self) -> &bam1_t { unsafe{self.inner.as_ref()} }
}

impl DerefMut for BamRec {
   #[inline]
   fn deref_mut(&mut self) -> &mut bam1_t {unsafe{ self.inner.as_mut() }}
}

impl AsRef<bam1_t> for BamRec {
   #[inline]
   fn as_ref(&self) -> &bam1_t { self}
}

impl AsMut<bam1_t> for BamRec {
   #[inline]
   fn as_mut(&mut self) -> &mut bam1_t { self}
}

unsafe impl Send for BamRec {}

impl HtsRead for BamRec {
   fn read(&mut self, hts: &mut Hts) -> io::Result<bool> {

      let (fp, hdr) = hts.hts_file_and_header();
      let hts_file = fp.as_mut();
      let res = if let Some(HtsHdr::Sam(hd)) = hdr {
         match unsafe { sam_read1(hts_file, hd.as_mut(), self.as_mut()) } {
            0..=c_int::MAX => Ok(true),
            -1 => Ok(false),
            _ => Err(hts_err("Error reading SAM/BAM record".to_string())),
         }
      } else {
         Err(hts_err(format!("Appropriate header not found for file {}", fp.name())))
      };
      res
   }

   fn read_itr(&mut self, hts: &mut Hts, itr: &mut HtsItr) -> io::Result<bool> {
      let fp = hts.hts_file();
      let i = if itr.multi() != 0 {
         unsafe { hts_itr_multi_next(fp.as_mut(), itr.as_mut(), self.as_mut() as *mut bam1_t as *mut c_void) }
      } else {
         let bgfp = if fp.as_ref().is_bgzf() != 0 {
            fp.as_ref().fp as *mut BGZF
         } else {
            null_mut::<BGZF>()
         };
         unsafe {hts_itr_next(bgfp, itr.as_mut(), self.as_mut() as *mut bam1_t as *mut c_void, fp.as_mut() as *mut htsFile as *mut c_void)}
      };
      if i >= 0 {
         Ok(true)
      } else if i == -1 {
         Ok(false)
      } else {
         Err(hts_err("Error reading record".to_string()))
      }
   }
}

impl Drop for BamRec {
   fn drop(&mut self) {
      unsafe { bam_destroy1(self.as_mut()) }
   }
}

pub(crate) fn check_tid(i: c_int) -> Option<usize> {
   if i >= 0 {
      Some(i as usize)
   } else {
      None
   }
}

impl BamRec {
   // Unlike the default behaviour in htslib, we initialize the new data field in
   // the bam1_t structure so we don't have to worry about this being zero
   pub fn new() -> io::Result<Self> {
      match NonNull::new(unsafe { bam_init1() }) {
         Some(mut b) => {
            let sz = 128;
            let data = unsafe { libc::calloc(sz, 1) as *mut c_char };
            if data.is_null() {
               Err(hts_err("Failed to allocate new BamRec".to_string()))
            } else {
               unsafe { b.as_mut().data = data }
               Ok(BamRec {
                  inner: b,
                  phantom: PhantomData,
               })
            }
         }
         None => Err(hts_err("Failed to allocate new BamRec".to_string())),
      }
   }

   pub fn copy(&self, bdst: &mut Self) -> io::Result<()> {
      match NonNull::new(unsafe { bam_copy1(bdst.as_mut(), self.as_ref()) }) {
         Some(b) => {
            bdst.inner = b;
            Ok(())
         }
         None => Err(hts_err("Failed to copy BamRec".to_string())),
      }
   }

   pub fn copy_from(b: Self) -> Self {
      BamRec {
         inner: b.inner,
         phantom: PhantomData,
      }
   }

   pub fn swap(&mut self, other: &mut Self) {
      std::mem::swap(&mut self.inner, &mut other.inner)
   }

   pub fn format(&mut self, hdr: &HtsHdr, s: &mut kstring_t) -> io::Result<()>{
      if let HtsHdr::Sam(h) = hdr {
         let e= self.as_mut().format(h, s);
         if e >= 0 {
            Ok(())
         } else {
            Err(hts_err(format!("Error formatting SAM record: err {}", e)))
         }
      } else {
         Err(hts_err("Wrong header type for SAM/BAM/CRAM format".to_string()))
      }
   }
}

fn aux_type2size(tp: u8) -> u8 {
   match tp {
      b'A' | b'c' | b'C' => 1,
      b's' | b'S' => 2,
      b'i' | b'I' | b'f' => 4,
      b'd' => 8,
      b'Z' | b'H' | b'B' => tp,
      _ => 0,
   }
}

pub struct BamAuxIter<'a> {
   data: &'a [u8],
}

impl<'a> Iterator for BamAuxIter<'a> {
   type Item = &'a [u8];
   fn next(&mut self) -> Option<Self::Item> {
      let ln = self.data.len();
      if ln < 3 {
         None
      } else {
         let mut l = 3;
         match aux_type2size(self.data[2]) {
            b'Z' | b'H' => {
               while l < ln && self.data[l] != 0 {
                  l += 1
               }
               if l < ln {
                  l += 1
               }
            }
            b'B' => {
               if ln - l < 5 {
                  return None;
               }
               let sz = aux_type2size(self.data[l]) as usize;
               let n =
                  u32::from_le_bytes(self.data[l + 1..l + 5].try_into().unwrap()) as usize;
               l += 5;
               if sz == 0 {
                  return None;
               }
               l += sz * n;
            }
            0 => return None,
            sz => l += sz as usize,
         }
         if l > ln {
            None
         } else {
            let (a, b) = self.data.split_at(l);
            self.data = b;
            Some(a)
         }
      }
   }
}

pub struct SeqQual(Box<[u8]>);

impl Deref for SeqQual {
   type Target = [u8];
   fn deref(&self) -> &[u8] {
      self.0.deref()
   }
}

impl DerefMut for SeqQual {
   fn deref_mut(&mut self) -> &mut [u8] {
      self.0.deref_mut()
   }
}

const FMT_BASES: [char; 256] = [
   'N', 'N', 'N', 'N', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G',
   'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C',
   'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A',
   'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
   'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G',
   'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C',
   'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A',
   'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
   'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G',
   'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C',
   'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A',
   'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
   'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G',
   'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
];

const FMT_QUAL: [u8; 256] = [
   33, 33, 33, 33, 34, 34, 34, 34, 35, 35, 35, 35, 36, 36, 36, 36, 37, 37, 37, 37, 38, 38, 38, 38,
   39, 39, 39, 39, 40, 40, 40, 40, 41, 41, 41, 41, 42, 42, 42, 42, 43, 43, 43, 43, 44, 44, 44, 44,
   45, 45, 45, 45, 46, 46, 46, 46, 47, 47, 47, 47, 48, 48, 48, 48, 49, 49, 49, 49, 50, 50, 50, 50,
   51, 51, 51, 51, 52, 52, 52, 52, 53, 53, 53, 53, 54, 54, 54, 54, 55, 55, 55, 55, 56, 56, 56, 56,
   57, 57, 57, 57, 58, 58, 58, 58, 59, 59, 59, 59, 60, 60, 60, 60, 61, 61, 61, 61, 62, 62, 62, 62,
   63, 63, 63, 63, 64, 64, 64, 64, 65, 65, 65, 65, 66, 66, 66, 66, 67, 67, 67, 67, 68, 68, 68, 68,
   69, 69, 69, 69, 70, 70, 70, 70, 71, 71, 71, 71, 72, 72, 72, 72, 73, 73, 73, 73, 74, 74, 74, 74,
   75, 75, 75, 75, 76, 76, 76, 76, 77, 77, 77, 77, 78, 78, 78, 78, 79, 79, 79, 79, 80, 80, 80, 80,
   81, 81, 81, 81, 82, 82, 82, 82, 83, 83, 83, 83, 84, 84, 84, 84, 85, 85, 85, 85, 86, 86, 86, 86,
   87, 87, 87, 87, 88, 88, 88, 88, 89, 89, 89, 89, 90, 90, 90, 90, 91, 91, 91, 91, 92, 92, 92, 92,
   93, 93, 93, 93, 94, 94, 94, 94, 95, 95, 95, 95, 96, 96, 96, 96,
];

impl fmt::Display for SeqQual {
   fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
      let mut s = String::with_capacity(self.len());
      if f.alternate() {
         for sq in self.iter() {
            s.push(FMT_QUAL[*sq as usize] as char)
         }
      } else {
         for sq in self.iter() {
            s.push(FMT_BASES[*sq as usize])
         }
      }
      f.write_str(&s)
   }
}

pub fn hts_sam_open_mode(mode: &mut [c_char], name: &str, fmt: Option<&str>) -> io::Result<()> {
   let res = if let Some(f) = fmt {
      unsafe {
         sam_open_mode(
            mode.as_mut_ptr(),
            get_cstr(name).as_ptr(),
            get_cstr(f).as_ptr(),
         )
      }
   } else {
      unsafe {
         sam_open_mode(
            mode.as_mut_ptr(),
            get_cstr(name).as_ptr(),
            std::ptr::null::<c_char>(),
         )
      }
   };

   if res != 0 {
      Err(hts_err(
         "sam_open_mode(): file mode not recognized".to_string(),
      ))
   } else {
      Ok(())
   }
}

