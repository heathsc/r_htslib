use std::{
   marker::PhantomData,
   ptr::{NonNull, null},
   ops::{Deref, DerefMut},
   ffi::{CString, CStr},
   io,
};

use libc::{c_char, c_int, c_void, size_t};

use crate::hts::HtsPos;
use super::{get_cstr, hts_err};

#[repr(C)]
pub struct regidx_t {
   _unused: [u8; 0],
}

#[repr(C)]
pub struct regitr_t {
   _unused: [u8; 0],
}

pub type RegidxParseF = unsafe extern fn(line: *const c_char, chr_beg: *mut *const c_char, chr_end: *mut *const c_char, beg: *mut HtsPos,
                                         end: *mut HtsPos, payload: *mut c_void, usr: *mut c_void) -> c_int;
pub type RegidxFreeF = unsafe extern fn(payload: *mut c_void);

#[link(name = "hts")]
extern "C" {
   fn regidx_destroy(idx: *mut regidx_t);
   fn regidx_init(fname: *const c_char, parse_f: Option<RegidxParseF>, free_f: Option<RegidxFreeF>, payload_size: size_t, usr: *mut c_void) -> *mut regidx_t;
   fn regidx_string(string: *const c_char, parse_f: Option<RegidxParseF>, free_f: Option<RegidxFreeF>, payload_size: size_t, usr: *mut c_void) -> *mut regidx_t;

   fn regidx_parse_bed(line: *const c_char, chr_beg: *mut *const c_char, chr_end: *mut *const c_char, beg: *mut HtsPos,
                       end: *mut HtsPos, payload: *mut c_void, usr: *mut c_void) -> c_int;
   fn regidx_parse_tab(line: *const c_char, chr_beg: *mut *const c_char, chr_end: *mut *const c_char, beg: *mut HtsPos,
                       end: *mut HtsPos, payload: *mut c_void, usr: *mut c_void) -> c_int;
   fn regidx_parse_reg(line: *const c_char, chr_beg: *mut *const c_char, chr_end: *mut *const c_char, beg: *mut HtsPos,
                       end: *mut HtsPos, payload: *mut c_void, usr: *mut c_void) -> c_int;
   fn regidx_parse_vcf(line: *const c_char, chr_beg: *mut *const c_char, chr_end: *mut *const c_char, beg: *mut HtsPos,
                       end: *mut HtsPos, payload: *mut c_void, usr: *mut c_void) -> c_int;

   fn regidx_seq_names(idx: *const regidx_t, n: *mut c_int) -> *const *const c_char;
   fn regidx_seq_nregs(idx: *const regidx_t, seq: *const c_char) -> c_int;
   fn regidx_nregs(idx: *const regidx_t) -> c_int;
   fn regidx_insert(idx: *mut regidx_t, line: *const c_char) -> c_int;
   fn regidx_push(idx: *mut regidx_t, chr_beg: *const c_char, chr_end: *const c_char, beg: HtsPos, end: HtsPos, payload: *const c_void) -> c_int;

}

pub struct RegIdx {
   inner: NonNull<regidx_t>,
   phantom: PhantomData<regidx_t>,
}

impl Deref for RegIdx {
   type Target = regidx_t;
   #[inline]
   fn deref(&self) -> &regidx_t { unsafe{self.inner.as_ref()} }
}

impl DerefMut for RegIdx {
   #[inline]
   fn deref_mut(&mut self) -> &mut regidx_t {unsafe{ self.inner.as_mut() }}
}

impl AsRef<regidx_t> for RegIdx {
   #[inline]
   fn as_ref(&self) -> &regidx_t { self}
}

impl AsMut<regidx_t> for RegIdx {
   #[inline]
   fn as_mut(&mut self) -> &mut regidx_t { self}
}

impl Drop for RegIdx {
   fn drop(&mut self) {
      unsafe { regidx_destroy(self.as_mut() ) }
   }
}

impl RegIdx {
   pub fn seq_names(&self) -> Vec<&CStr> {
      let mut n: c_int = 0;
      let p = unsafe { regidx_seq_names(self.as_ref(), &mut n as *mut c_int) };

      // Check for validity
      assert!(n == 0 || (n > 0 && !p.is_null()));

      // Build vector of CStr
      let mut v = Vec::with_capacity(n as usize);
      if n > 0 {
         let names = unsafe {std::slice::from_raw_parts(p, n as usize)};
         for s in names {
            v.push(unsafe { CStr::from_ptr(*s) })
         }
      }
      v
   }

   pub fn seq_nregs<S: AsRef<str>>(&self, s: S) -> c_int {
      unsafe{ regidx_seq_nregs(self.as_ref(), get_cstr(s).as_ptr()) }
   }

   pub fn nregs(&self) -> c_int {
      unsafe{ regidx_nregs(self.as_ref()) }
   }

   pub fn insert<S: AsRef<str>>(&mut self, line: S) -> io::Result<()> {
      let line = get_cstr(line);
      if unsafe { regidx_insert(self.as_mut(), line.as_ptr()) } == 0 {
         Ok(())
      } else {
         Err(hts_err(format!("Error inserting region {:?}", line)))
      }
   }

   pub fn insert_list<S: AsRef<str>>(&mut self, line: S, delim: char) -> io::Result<()> {
      for s in line.as_ref().split(delim) {
         self.insert(s)?
      }
      Ok(())
   }

   pub fn insert_push<S: AsRef<str>>(&mut self, chr: S, beg: Option<usize>, end: Option<usize>, payload: Option<*mut c_void>) -> io::Result<()> {
      let chr = chr.as_ref();
      let chr_len = chr.len();
      if chr_len == 0 { return Err(hts_err("RegIdx::insert_push() called with empty chromosome string".to_string()))}
      let beg = beg.map(|x| x as HtsPos).unwrap_or(0);
      let end = end.map(|x| x as HtsPos).unwrap_or(HtsPos::MAX);
      let chr_beg = chr.as_ptr() as *const c_char;
      let chr_end = unsafe{ chr_beg.add(chr_len - 1) };
      let payload = payload.unwrap_or(null::<*mut c_void>() as *mut c_void);
      if unsafe { regidx_push(self.as_mut(), chr_beg, chr_end, beg, end, payload)} == 0 {
         Ok(())
      } else {
         Err(hts_err(format!("Error pushing region {}", chr)))
      }
   }

   pub fn insert_chr<S: AsRef<str>>(&mut self, chr: S) -> io::Result<()> {
      self.insert_push(chr, None, None, None)
   }
}

#[derive(Debug)]
pub enum RegIdxInitInput {
   Fname(CString),
   String(CString),
   None,
}

impl Default for RegIdxInitInput {
   fn default() -> Self { Self::None }
}

#[derive(Debug)]
pub enum ParseF {
   ParseBed,
   ParseTab,
   ParseReg,
   ParseVcf,
   Custom(RegidxParseF),
   None,
}

impl Default for ParseF {
   fn default() -> Self { Self::None }
}

impl ParseF {
   pub fn parse_f(&self) -> Option<RegidxParseF> {
      match self {
         Self::ParseBed => Some(regidx_parse_bed),
         Self::ParseTab => Some(regidx_parse_tab),
         Self::ParseReg => Some(regidx_parse_reg),
         Self::ParseVcf => Some(regidx_parse_vcf),
         Self::Custom(f) => Some(*f),
         Self::None => None,
      }
   }
}

#[derive(Default, Debug)]
pub struct RegIdxInit {
   input: RegIdxInitInput,
   parse_f: ParseF,
   free_f: Option<RegidxFreeF>,
   payload_size: size_t,
   usr: Option<*mut c_void>,
}

impl RegIdxInit {
   pub fn new() -> Self { Self::default() }

   pub fn fname<S: AsRef<str>>(mut self, s: S) -> Self {
      self.input = RegIdxInitInput::Fname(get_cstr(s));
      self
   }

   pub fn string<S: AsRef<str>>(mut self, s: S) -> Self {
      self.input = RegIdxInitInput::String(get_cstr(s));
      self
   }

   pub fn parse_f(mut self, parse_f: ParseF) -> Self {
      self.parse_f = parse_f;
      self
   }

   pub fn build(self) -> Option<RegIdx> {

      let usr= self.usr.unwrap_or(null::<*mut c_void>() as *mut c_void);
      NonNull::new( match self.input {
         RegIdxInitInput::Fname(s) => unsafe { regidx_init(s.as_ptr(), self.parse_f.parse_f(), self.free_f, self.payload_size, usr) },
         RegIdxInitInput::String(s) => unsafe { regidx_string(s.as_ptr(), self.parse_f.parse_f(), self.free_f, self.payload_size, usr) },
         RegIdxInitInput::None => unsafe { regidx_string(null(), self.parse_f.parse_f(), self.free_f, self.payload_size, usr) },
      }).map(|inner| RegIdx{inner, phantom: PhantomData})
   }
}