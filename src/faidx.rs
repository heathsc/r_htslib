use std::{
    io,
    ptr::{NonNull, null},
    marker::PhantomData,
    ops::{Deref, DerefMut},
    path::Path,
    os::unix::ffi::OsStrExt,
    ffi::CString,
};

use super::{from_cstr, get_cstr, hts_err, HtsPos};
use libc::{c_char, c_int, c_void, free};

#[repr(C)]
pub struct faidx_t {
    _unused: [u8; 0],
}

impl faidx_t {
    pub fn nseq(&self) -> usize {
        let l = unsafe{ faidx_nseq(self)};
        l as usize
    }
    pub fn iseq(&self, i: usize) -> &str {
        if i > self.nseq() { panic!("Sequence ID {} out of range", i); }
        from_cstr(unsafe { faidx_iseq(self, i as libc::c_int) })
    }

    pub fn seq_len<S: AsRef<str>>(&self, cname: S) -> Option<usize> {
        let cname = cname.as_ref();
        let len = unsafe{ faidx_seq_len(self, get_cstr(cname).as_ptr())};
        if len < 0 { None } else { Some(len as usize) }
    }

    // Attempts to load reference sequence from file
    // x and y are 1 offset coordinates.  Setting x to 0 or 1 will load from the start of the contig.  Setting y to None
    // or a very large value will load until the end of the chromosome.
    // Returns errors if the chromosome is not found, the coordinates are invalid (i.e., y < x) or an IO error occurred
    pub fn fetch_seq<S: AsRef<str>>(&self, cname: S, x: usize, y: Option<usize>) -> io::Result<Sequence> {
        let cname = cname.as_ref();
        if let Some(seq_len) = self.seq_len(cname) {
            let y = y.map(|z| z.min(seq_len)).unwrap_or(seq_len);
            let x = x.saturating_sub(1);
            if y <= x { Err(hts_err(format!("Requested Sequence {}:{}-{} has zero or negative length", cname, x, y))) }
            else {
                let mut len: HtsPos = 0;
                let seq = unsafe{ faidx_fetch_seq64(self, get_cstr(cname).as_ptr(), x as HtsPos, (y - 1) as HtsPos, &mut len) };
                if len == -2 { Err(hts_err(format!("Sequence {} not found", cname))) }
                else if len < 0 || seq.is_null() { Err(hts_err(format!("Loading of sequence data for {} failed", cname))) }
                else { Ok(Sequence{inner: NonNull::new(seq as *mut u8).unwrap(), phantom: PhantomData, start: x + 1, len: len as usize}) }
            }
        } else { Err(hts_err(format!("Sequence {} not found", cname))) }
    }
}

extern "C" {
    fn fai_load(fn_: *const c_char) -> *mut faidx_t;
    fn fai_load3(fn_: *const c_char, fnai: *const c_char, fngzi: *const c_char, flags: c_int) -> *mut faidx_t;
    fn faidx_nseq(fai: *const faidx_t) -> c_int;
    fn faidx_iseq(fai: *const faidx_t, n: c_int) -> *const c_char;
    fn faidx_seq_len(fai: *const faidx_t, seq: *const c_char) -> c_int;
    fn faidx_fetch_seq64(
        fai: *const faidx_t,
        cname: *const c_char,
        x: HtsPos,
        y: HtsPos,
        len: *mut HtsPos,
    ) -> *mut c_char;
}

pub struct Faidx {
    inner: NonNull<faidx_t>,
    phantom: PhantomData<faidx_t>,
}

unsafe impl Send for Faidx {}
unsafe impl Sync for Faidx {}

impl Deref for Faidx {
    type Target = faidx_t;
    #[inline]
    fn deref(&self) -> &faidx_t { unsafe{self.inner.as_ref()} }
}

impl DerefMut for Faidx {
    #[inline]
    fn deref_mut(&mut self) -> &mut faidx_t {unsafe{ self.inner.as_mut() }}
}

impl Faidx {
    pub fn load<S: AsRef<Path>>(name: S) -> io::Result<Faidx> {
        let cname = CString::new(name.as_ref().as_os_str().as_bytes())?;

        match NonNull::new(unsafe{ fai_load3(cname.as_ptr(), null(), null(), 0)}) {
            None => Err(hts_err(format!("Failed to load reference file index {}", name.as_ref().display()))),
            Some(idx) => Ok(Faidx{inner: idx, phantom: PhantomData}),
        }
    }
    pub fn load_or_create<S: AsRef<Path>>(name: S) -> io::Result<Faidx> {
        let cname = CString::new(name.as_ref().as_os_str().as_bytes())?;

        match NonNull::new(unsafe{ fai_load(cname.as_ptr())}) {
            None => Err(hts_err(format!("Failed to load reference file index {}", name.as_ref().display()))),
            Some(idx) => Ok(Faidx{inner: idx, phantom: PhantomData}),
        }
    }
}

pub struct Sequence {
    inner: NonNull<u8>,
    phantom: PhantomData<u8>,
    start: usize,
    len: usize,
}

unsafe impl Send for Sequence {}
unsafe impl Sync for Sequence {}

impl Drop for Sequence {
    fn drop(&mut self) {
        unsafe { free(self.inner.as_ptr() as *mut c_void) }
    }
}

impl Sequence {
    // Get sequence between x and y inclusive (1 offset)
    pub fn get_seq(&self, x: usize, y: usize) -> io::Result<&[u8]> {
        if x < 1 || x < self.start || x > y {
            Err(hts_err("Invalid coordinates".to_string()))
        } else {
            let a = x - self.start;
            let b = (y + 1 - self.start).min(self.len);
            let slice = self.seq();
            Ok(&slice[a..b])
        }
    }

    // Get entire loaded sequence as a slice
    pub fn seq(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.inner.as_ptr(), self.len) }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

