use std::ffi::{CStr, CString};
use std::io;
use std::io::{Error, ErrorKind};

pub mod hts;
pub use hts::*;
pub mod sam;
pub use sam::*;
pub mod vcf;
pub use vcf::*;
pub mod faidx;
pub use faidx::*;
pub mod kstring;
pub mod tbx;
pub use tbx::*;
pub mod regidx;
pub use kstring::*;

#[inline]
pub fn get_cstr<S: AsRef<str>>(s: S) -> CString {
    CString::new(s.as_ref().as_bytes()).unwrap()
}

#[inline]
pub fn try_from_cstr<'a>(c: *const i8) -> Option<&'a str> {
    if c.is_null() {
        None
    } else {
        Some(from_cstr(c))
    }
}

#[inline]
pub fn from_cstr<'a>(c: *const i8) -> &'a str {
    c_to_cstr(c)
        .to_str()
        .expect("String not UTF8")
}

#[inline]
pub fn c_to_cstr<'a>(c: *const i8) -> &'a CStr {
    if c.is_null() {
        panic!("from_cstr() called with a NULL");
    }
    unsafe {
        CStr::from_ptr(c)
    }
}

pub fn hts_err(s: String) -> io::Error {
    Error::new(ErrorKind::Other, s)
}
