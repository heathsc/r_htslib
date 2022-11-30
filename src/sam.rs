use std::convert::TryInto;
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};
use std::ptr::{NonNull, null_mut};

use std::{fmt, io};

use crate::{
   get_cstr, htsFile, hts_err, kstring_t, Hts, HtsFile,
   HtsHdr, HtsItr, HtsRead, HtsWrite, BGZF, hts_itr_next, hts_itr_multi_next,
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

unsafe impl Send for SamHeader {}
unsafe impl Sync for SamHeader {}

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

   pub fn seq_names(&self) -> Vec<&str> {
      let n = self.nref();
      (0..n).map(|i| self.tid2name(i)).collect()
   }

   pub fn seq_lengths(&self) -> Vec<usize> {
      let n = self.nref();
      (0..n).map(|i| self.tid2len(i)).collect()
   }
}

impl Clone for SamHeader {
   fn clone(&self) -> Self {
      self.dup().expect("Error duplicating SAM/BAM header")
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
unsafe impl Sync for BamRec {}

impl HtsRead for BamRec {
   fn read(&mut self, hts: &mut Hts) -> io::Result<bool> {

      let (fp, hdr) = hts.hts_file_and_header();
      let hts_file = fp.as_mut();
      let res = if let Some(HtsHdr::Sam(hd)) = hdr {
         match unsafe { sam_read1(hts_file, hd.as_mut(), self.as_mut()) } {
            0..=c_int::MAX => Ok(true),
            -1 => Ok(false),
            _ => Err(hts_err("Error reading SAM/BAM/CRAM record".to_string())),
         }
      } else {
         Err(hts_err(format!("Appropriate header not found for file {}", fp.name())))
      };
      res
   }

   fn read_itr(&mut self, hts: &mut Hts, itr: &mut HtsItr) -> io::Result<bool> {
      let fp = hts.hts_file_mut();
      let i = if itr.multi() != 0 {
         unsafe { hts_itr_multi_next(fp.as_mut(), itr.as_mut(), self.as_mut() as *mut bam1_t as *mut c_void) }
      } else {
         let bgfp = if fp.as_ref().is_bgzf() != 0 {
            unsafe {fp.as_ref().fp.bgzf}
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

impl HtsWrite for BamRec {
   fn write(&mut self, hts: &mut Hts) -> io::Result<()> {
      let (fp, hdr) = hts.hts_file_and_header();
      let hts_file = fp.as_mut();
      let res = if let Some(HtsHdr::Sam(hd)) = hdr {
         match unsafe { sam_write1(hts_file, hd.as_mut(), self.as_mut()) } {
            0..=c_int::MAX => Ok(()),
            _ => Err(hts_err("Error writing SAM/BAM/CRAM record".to_string())),
         }
      } else {
         Err(hts_err(format!("Appropriate header not found for file {}", fp.name())))
      };
      res
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

enum AuxType {
   Size(usize),
   NullTerm,
   Array,
   None,
}

impl AuxType {
   fn from_u8_code(tp: u8) -> Self {
      match tp {
         b'A' | b'c' | b'C' => Self::Size(1),
         b's' | b'S' => Self::Size(2),
         b'i' | b'I' | b'f' => Self::Size(4),
         b'd' => Self::Size(8),
         b'Z' | b'H' => Self::NullTerm,
         b'B' => Self::Array,
         _ => Self::None,
      }
   }
}

/*
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
*/

// Returns length of c-type null terminated u8 slice including the null
// If no null is found, returns size of input slice
fn null_term_length(p: &[u8]) -> usize {
   for (i, c) in p.iter().enumerate() {
      if *c == 0 {
         return i + 1
      }
   }
   p.len()
}

// Returns length of Bam aux tag value (i.e., not include the 2 character tag or the initial one character type)
fn tag_length(p: &[u8]) -> Option<usize> {
   let l = p.len();
   if l < 2 {
      None
   } else {
      let chk = |x| if x >= l { None } else { Some(x) };

      match AuxType::from_u8_code(p[0]) {
         AuxType::Size(x) => chk(x),
         AuxType::NullTerm => Some(null_term_length(&p[1..])),
         AuxType::Array => {
            if l < 6 { None } else {
               match AuxType::from_u8_code(p[1]) {
                  AuxType::Size(x) => {
                     let n =
                         u32::from_le_bytes(p[2..6].try_into().unwrap()) as usize;
                     chk(x * n)
                  },
                  _ => None,
               }
            }
         }
         AuxType::None => None,
      }
   }
}

pub trait BamAux<'a> {
   fn data(&self) -> &[u8];

   fn tag(&self) -> &[u8] {
      &self.data()[..2]
   }
   fn val_type(&self) -> u8 {
      self.data()[2]
   }
   fn value(&self) -> &[u8] {
      &self.data()[3..]
   }
   fn len(&self) -> usize { self.data().len() }

   fn is_empty(&self) -> bool { self.data().is_empty() }

   fn as_ptr(&self) -> *const u8 {
      self.data().as_ptr()
   }
}

pub struct BamAuxItemMut<'a> {
   data: &'a mut [u8],
}

impl <'a> BamAux<'a> for BamAuxItemMut<'a> {
   fn data(&self) -> &[u8] {
      self.data
   }
}

impl <'a> BamAuxItemMut <'a> {
   pub fn as_mut_ptr(&mut self) -> *mut u8 {
      self.data.as_mut_ptr()
   }
}

pub struct BamAuxIterMut<'a> {
   ptr: NonNull<u8>,
   size: usize,
   _marker: PhantomData<&'a mut u8>
}

impl <'a> BamAuxIterMut<'a> {
   fn new(p: *mut u8, size: usize) -> Option<Self> {
      NonNull::new(p).map(|ptr| Self { ptr, size, _marker: PhantomData })
   }
}

impl<'a> Iterator for BamAuxIterMut<'a> {
   type Item = BamAuxItemMut<'a>;
   fn next(&mut self) -> Option<Self::Item> {
      let ln = self.size;
      if ln < 3 {
         None
      } else {
         let data = unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.size) };
         tag_length(&data[2..]).map(|x| {
            let p1 = self.ptr.as_ptr();
            let l1 = x + 3;
            let p2 = unsafe { p1.add(x + 3) };
            assert!(ln >= x + 3);
            self.size -= x + 3;
            self.ptr = NonNull::new(p2).unwrap();
            let d = unsafe { std::slice::from_raw_parts_mut(p1, l1) };
            BamAuxItemMut{data: d}
         })
      }
   }
}

pub struct BamAuxIter<'a> {
   data: &'a [u8],
}

impl <'a> BamAuxIter<'a> {
   pub fn new(data: &'a [u8]) -> Self {
      Self {data}
   }
}

impl <'a> BamAux<'a> for BamAuxItem<'a> {
   fn data(&self) -> &[u8] {
      self.data
   }
}

pub struct BamAuxItem<'a> {
   data: &'a [u8],
}

impl<'a> Iterator for BamAuxIter<'a> {
   type Item = BamAuxItem<'a>;
   fn next(&mut self) -> Option<Self::Item> {
      let ln = self.data.len();
      if ln < 3 {
         None
      } else {
         tag_length(&self.data[2..]).map(|x| {
            let (a, b) = self.data.split_at(x + 3);
            self.data = b;
            BamAuxItem{data: a}
         })
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

