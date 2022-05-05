pub mod lib;
pub use lib::*;

use libc::{c_void, c_int};

use std::{
    io,
    marker::PhantomData,
    ptr::{null, NonNull},
    ops::{Deref, DerefMut},
    ffi::CStr,
};

use crate::{
    tbx_index_load3, VcfHeader, SamHeader, Tbx, vcf_read_itr, bcf_read_itr, tbx_read_itr, tbx_name2id,
    sam_itr_queryi, sam_hdr_t, bam_name2id, BamRec, BcfRec, TbxRec, bcf_hdr_name2id,
    get_cstr, hts_err, bcf_hdr_t, tbx_t
};

pub struct Hts {
    pub(crate) hts_file: HtsFile,
    pub(crate) header: Option<HtsHdr>,
    pub(crate) tbx: Option<Tbx>,
    pub(crate) idx: Option<HtsIndex>,
}

unsafe impl Send for Hts {}

impl Hts {
    pub fn open<S: AsRef<str>>(name: S, mode: &str) -> io::Result<Self> {
        Self::open_format_(name, mode, None)
    }

    pub fn open_format<S: AsRef<str>>(name: S, mode: &str, fmt: &HtsFormat) -> io::Result<Self> {
        Self::open_format_(name, mode, Some(fmt))
    }

    fn open_format_<S: AsRef<str>>(name: S, mode: &str, fmt: Option<&HtsFormat>) -> io::Result<Self> {
        let name = name.as_ref();
        let mut fp = HtsFile::open_format_(name, mode, fmt)?;
        let (header, tbx) = if fp.as_ref().is_write() == 0 {
            match fp.format().format {
                htsExactFormat::Sam | htsExactFormat::Bam | htsExactFormat::Cram => (Some(HtsHdr::Sam(SamHeader::read(&mut fp)?)), None),
                htsExactFormat::Bcf | htsExactFormat::Vcf => (Some(HtsHdr::Vcf(VcfHeader::read(&mut fp)?)), None),
                _ => if let Ok(tbx) = Tbx::load(name) { (None, Some(tbx)) } else { (None, None) }
            }
        } else { (None, None) };

        Ok(Self {
            hts_file: fp,
            header,
            tbx,
            idx: None,
        })
    }

    pub fn header(&self) -> Option<&HtsHdr> { self.header.as_ref() }

    pub fn header_mut(&mut self) -> Option<&mut HtsHdr> { self.header.as_mut() }

    pub fn tbx(&self) -> Option<&Tbx> { self.tbx.as_ref() }

    pub fn tbx_mut(&mut self) -> Option<&mut Tbx> { self.tbx.as_mut() }

    pub fn hts_file(&mut self) -> &mut HtsFile { &mut self.hts_file }

    pub fn hts_file_and_header(&mut self) -> (&mut HtsFile, Option<&mut HtsHdr>) {
        (&mut self.hts_file, self.header.as_mut())
    }

    pub fn hts_file_and_tbx(&mut self) -> (&mut HtsFile, Option<&mut Tbx>) {
        (&mut self.hts_file, self.tbx.as_mut())
    }

    pub fn index_load(&mut self) -> io::Result<()> {
        self.index_load3(None, 0)
    }

    pub fn index_load2(&mut self, fnidx: Option<&str>) -> io::Result<()> {
        self.index_load3(fnidx, 0)
    }

    pub fn name(&self) -> &str { self.hts_file.name() }

    pub fn has_index(&self) -> bool {
        self.idx.is_some() || self.tbx.is_some()
    }

    pub fn index_load3(&mut self, fnidx: Option<&str>, flags: usize) -> io::Result<()> {
        let flags = (flags as c_int) | HTS_IDX_SILENT_FAIL;
        let fname = self.hts_file.name_ptr();
        let fnidx = fnidx.map(get_cstr);
        let have_fnidx = fnidx.is_some();
        let fnidx_ptr = match fnidx {
            Some(s) => s.as_ptr(),
            None => null(),
        };

        let mk_hts_idx = |p: *mut hts_idx_t| {
            NonNull::new(p).map(|q| HtsIndex { inner: q, phantom: PhantomData})
        };

        let idx = match self.hts_file.format().format {
            htsExactFormat::Bam | htsExactFormat::Cram => {
                self.idx = mk_hts_idx(unsafe { sam_index_load3(self.hts_file.as_mut(), fname, fnidx_ptr, flags) });
                self.idx.is_some()
            },
            htsExactFormat::Sam => {
                self.idx = mk_hts_idx(unsafe { sam_index_load3(self.hts_file.as_mut(), fname, fnidx_ptr, flags) });
                self.idx.is_some() || {
                    if let Some(q) = NonNull::new(unsafe { tbx_index_load3(fname, fnidx_ptr, flags) }) {
                        self.tbx = Some(Tbx::new(q));
                        true
                    } else { false }
                }
            },
            htsExactFormat::Bcf => {
                self.idx = mk_hts_idx(unsafe { hts_idx_load3(fname, fnidx_ptr, HTS_FMT_CSI, flags)});
                self.idx.is_some()
            },
            htsExactFormat::Vcf => {
                self.idx = mk_hts_idx(unsafe { hts_idx_load3(fname, fnidx_ptr, HTS_FMT_CSI, flags)});
                self.idx.is_some() || {
                    if let Some(q) = NonNull::new(unsafe { tbx_index_load3(fname, fnidx_ptr, flags) }) {
                        self.tbx = Some(Tbx::new(q));
                        true
                    } else { false }
                }
            },
            _ => {
                // Retry to load index for Tabix file if we have been supplied with an index name
                if have_fnidx {
                    if let Some(q) = NonNull::new(unsafe { tbx_index_load3(fname, fnidx_ptr, flags) }) {
                        self.tbx = Some(Tbx::new(q));
                    }
                }
                self.tbx.is_some()
            }
        };
        if idx { Ok(()) } else { Err(hts_err(format!("Failed to load index for file {}", self.name()))) }
    }

    pub fn index(&self) -> Option<&hts_idx_t> {
        if self.has_index() {
            self.idx.as_ref().map_or_else(
                || self.tbx.as_ref().map(|tbx| tbx.idx()),
                |p| Some(p.as_ref()))
        } else { None }
    }

    pub fn rec_type(&self) -> Option<HtsRecType> {
        if self.tbx.is_some() { Some(HtsRecType::Tbx) } else {
            self.header().map(|h| match h {
                HtsHdr::Sam(_) => HtsRecType::Sam,
                HtsHdr::Vcf(_) => HtsRecType::Vcf,
            })
        }
    }

    pub fn itr_queryi(&mut self, tid: usize, begin: HtsPos, end: HtsPos) -> io::Result<HtsItr> {
        if !self.has_index() { self.index_load()? }
        let idx= self.index().unwrap();
        let hdr = self.header.as_ref();

        if let Some(itr) = NonNull::new(match (hdr, self.tbx.as_ref()) {
            (_, Some(_)) => unsafe { hts_itr_query(idx, tid as c_int, begin, end, tbx_read_itr) },
            (Some(HtsHdr::Sam(_)), _) => unsafe { sam_itr_queryi(idx, tid as c_int, begin, end) },
            (Some(HtsHdr::Vcf(_)), _) if matches!(self.hts_file.format().format, htsExactFormat::Bcf) => unsafe { hts_itr_query(idx, tid as c_int, begin, end, bcf_read_itr) },
            (Some(HtsHdr::Vcf(_)), _) => unsafe { hts_itr_query(idx, tid as c_int, begin, end, vcf_read_itr) },
            (None, None) => return Err(hts_err(format!("Error making iterator for file {}", self.name()))),
        }) {
            Ok(HtsItr { inner: itr, phantom: PhantomData })
        } else {
            Err(hts_err(format!("Error making iterator for file {}", self.name())))
        }
    }

    pub fn itr_querys(&mut self, reg: &CStr) -> io::Result<HtsItr> {
        match reg.to_bytes() {
            [b'.'] => self.itr_queryi(HTS_IDX_START as usize, 0, 0),
            [b'*'] => self.itr_queryi(HTS_IDX_NOCOOR as usize, 0, 0),
            _ => {
                let (tid, begin, end) = self.parse_region(reg)?;
                self.itr_queryi(tid as usize, begin, end)
            },
        }
    }

    fn parse_region(&mut self, reg: &CStr) -> io::Result<(c_int, HtsPos, HtsPos)> {
        let mut beg: HtsPos = 0;
        let mut end: HtsPos = 0;
        let mut tid: c_int = 0;

        let (get_id, hdr):(HtsName2Id, *mut c_void) = if let Some(tbx) = self.tbx.as_mut() {
            // Tabix file
            (tbx_name2id, tbx.as_mut() as *mut tbx_t as *mut c_void)
        } else {
            // VCF/BCF or SAM/BAM/CRAM files
            match self.header_mut() {
                Some(HtsHdr::Vcf(h)) => (bcf_hdr_name2id, h.as_mut() as *mut bcf_hdr_t as *mut c_void),
                Some(HtsHdr::Sam(h)) => (bam_name2id, h.as_mut() as *mut sam_hdr_t as *mut c_void),
                _ => return Err(hts_err(format!("Error making iterator for file {}", self.name()))),
            }
        };

        if unsafe {
            hts_parse_region(reg.as_ptr(), &mut tid, &mut beg, &mut end, get_id, hdr, HTS_PARSE_THOUSANDS_SEP)
        }.is_null() { return Err(hts_err(format!("Error parsing region {:?}", reg))) }

        Ok((tid, beg, end))
    }

    pub fn reader<T>(&mut self) -> HtsReader<T> {
        HtsReader {
            hts: self,
            _phantom: PhantomData,
        }
    }

    pub fn itr_reader<T>(&mut self, itr: HtsItr) -> HtsItrReader<T> {
        HtsItrReader {
            hts: self,
            itr,
            _phantom: PhantomData,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum HtsRecType {
    Sam, Vcf, Tbx
}

pub enum HtsRec {
    Sam(BamRec),
    Vcf(BcfRec),
    Tbx(TbxRec),
}

pub struct HtsFile {
    pub(crate) inner: NonNull<htsFile>,
    phantom: PhantomData<htsFile>,
}

unsafe impl Send for HtsFile {}

impl Deref for HtsFile {
    type Target = htsFile;
    #[inline]
    fn deref(&self) -> &htsFile { unsafe{self.inner.as_ref()} }
}

impl DerefMut for HtsFile {
    #[inline]
    fn deref_mut(&mut self) -> &mut htsFile {unsafe{ self.inner.as_mut() }}
}

impl AsRef<htsFile> for HtsFile {
    #[inline]
    fn as_ref(&self) -> &htsFile { self}
}

impl AsMut<htsFile> for HtsFile {
    #[inline]
    fn as_mut(&mut self) -> &mut htsFile { self}
}

impl Drop for HtsFile {
    fn drop(&mut self) {
        unsafe { hts_close(self.as_mut()) };
    }
}

impl HtsFile {
    pub fn open<S: AsRef<str>>(name: S, mode: &str) -> io::Result<Self> {
        Self::open_format_(name, mode, None)
    }

    pub fn open_format<S: AsRef<str>>(name: S, mode: &str, fmt: &HtsFormat) -> io::Result<Self> {
        Self::open_format_(name, mode, Some(fmt))
    }

    fn open_format_<S: AsRef<str>>(name: S, mode: &str, fmt: Option<&HtsFormat>) -> io::Result<Self> {
        let name = name.as_ref();
        let mode = get_cstr(mode);
        let fmt = match fmt {
            Some(f) => f.inner.as_ptr(),
            None => null(),
        };
        if let Some(hfile) = NonNull::new(unsafe { hts_open_format(get_cstr(name).as_ptr(), mode.as_ptr(), fmt) }) {
            assert!(!unsafe {hfile.as_ref()}.fn_.is_null());
            let hfile = hfile.cast::<htsFile>();
            Ok(HtsFile {
                inner: hfile,
                phantom: PhantomData,
            })
        } else {
            Err(hts_err(format!("Couldn't open file {}", name)))
        }
    }

    pub fn set_thread_pool(&mut self, tp: &HtsThreadPool) -> io::Result<()> {
        if unsafe{hts_set_thread_pool(self.as_mut(), tp)} == 0 {
            Ok(())
        } else {
            Err(hts_err(format!("Failed to set thread pool for file {}", self.name())))
        }
    }
}

#[derive(Debug)]
pub enum HtsHdr {
    Vcf(VcfHeader),
    Sam(SamHeader),
}

#[repr(C)]
pub struct HtsThreadPool {
    inner: NonNull<hts_tpool>,
    qsize: c_int,
    phantom: PhantomData<hts_tpool>,
}

unsafe impl Send for HtsThreadPool {}

impl Drop for HtsThreadPool {
    fn drop(&mut self) {
        unsafe { hts_tpool_destroy(self.inner.as_ptr() ) }
    }
}

impl HtsThreadPool {
    pub fn new(nthreads: usize) -> Option<Self> {
        NonNull::new(unsafe {hts_tpool_init(nthreads as c_int)})
           .map(|tpool| Self {
                inner: tpool,
                qsize: 0,
                phantom: PhantomData,
            })
    }
}

pub struct HtsIndex {
    inner: NonNull<hts_idx_t>,
    phantom: PhantomData<hts_idx_t>,
}

unsafe impl Send for HtsIndex {}

impl Deref for HtsIndex {
    type Target = hts_idx_t;
    #[inline]
    fn deref(&self) -> &hts_idx_t { unsafe{self.inner.as_ref()} }
}

impl DerefMut for HtsIndex {
    #[inline]
    fn deref_mut(&mut self) -> &mut hts_idx_t {unsafe{ self.inner.as_mut() }}
}

impl AsRef<hts_idx_t> for HtsIndex {
    #[inline]
    fn as_ref(&self) -> &hts_idx_t { self}
}

impl AsMut<hts_idx_t> for HtsIndex {
    #[inline]
    fn as_mut(&mut self) -> &mut hts_idx_t { self}
}

impl Drop for HtsIndex {
    fn drop(&mut self) {
        unsafe{ hts_idx_destroy(self.as_mut()) }
    }
}

pub struct HtsFormat {
    inner: NonNull<htsFormat>,
    phantom: PhantomData<htsFormat>,
    owned: bool,
}

impl Drop for HtsFormat {
    fn drop(&mut self) {
        if self.owned {
            let fmt = unsafe { Box::from_raw(self.inner.as_ptr()) };
            drop(fmt)
        }
    }
}

impl Default for HtsFormat {
    fn default() -> Self {
        let fmt = Box::new(htsFormat::default());
        let fmt: *mut htsFormat = Box::into_raw(fmt);
        Self {
            inner: NonNull::new(fmt).unwrap(),
            phantom: PhantomData,
            owned: true,
        }
    }
}

impl HtsFormat {
    pub fn parse_format<S: AsRef<str>>(s: S) -> io::Result<Self> {
        let s = s.as_ref();
        let fmt = Box::new(htsFormat::default());
        let fmt: *mut htsFormat = Box::into_raw(fmt);
        if unsafe { hts_parse_format(fmt, get_cstr(s).as_ptr()) } == 0 {
            Ok(Self {
                inner: NonNull::new(fmt).unwrap(),
                phantom: PhantomData,
                owned: true,
            })
        } else {
            Err(hts_err(format!("Error parsing '{}' as hts format", s)))
        }
    }

    pub fn opt_add<S: AsRef<str>>(&mut self, s: S) -> io::Result<()> {
        let s = s.as_ref();
        if unsafe { hts_opt_add(&mut self.inner.as_mut().specific, get_cstr(s).as_ptr()) } != 0 {
            Err(hts_err(format!("hts_add_opt(): Illegal option {}", s)))
        } else {
            Ok(())
        }
    }

    pub fn is_compressed(&self) -> bool {
        let format = unsafe { &self.inner.as_ref() };
        format.compression == htsCompression::Bgzf || format.format == htsExactFormat::Cram
    }
}

pub struct HtsItr{
    inner: NonNull<hts_itr_t>,
    phantom: PhantomData<hts_itr_t>,
}

impl Deref for HtsItr {
    type Target = hts_itr_t;
    #[inline]
    fn deref(&self) -> &hts_itr_t { unsafe{self.inner.as_ref()} }
}

impl DerefMut for HtsItr {
    #[inline]
    fn deref_mut(&mut self) -> &mut hts_itr_t {unsafe{ self.inner.as_mut() }}
}

impl AsRef<hts_itr_t> for HtsItr {
    #[inline]
    fn as_ref(&self) -> &hts_itr_t { self}
}

impl AsMut<hts_itr_t> for HtsItr {
    #[inline]
    fn as_mut(&mut self) -> &mut hts_itr_t { self}
}

impl Drop for HtsItr {
    fn drop(&mut self) {
        unsafe { hts_itr_destroy(self.as_mut()) };
    }
}

pub trait HtsRead {
    fn read(&mut self, hts: &mut Hts) -> io::Result<bool>;
    fn read_itr(&mut self, hts: &mut Hts, itr: &mut HtsItr) -> io::Result<bool>;
}

pub struct HtsReader<'a, T> {
    hts: &'a mut Hts,
    _phantom: PhantomData<T>,
}

impl <'a, T: HtsRead> HtsReader<'a, T> {
    pub fn read(&mut self, rec: &mut T) -> io::Result<bool> {
        rec.read(self.hts)
    }

    pub fn header(&self) -> Option<&HtsHdr> {
        self.hts.header()
    }
}

pub struct HtsItrReader<'a, T> {
    hts: &'a mut Hts,
    itr: HtsItr,
    _phantom: PhantomData<T>,
}

impl <'a, T: HtsRead> HtsItrReader<'a, T> {
    pub fn read(&mut self, rec: &mut T) -> io::Result<bool> {
        rec.read_itr(self.hts, &mut self.itr)
    }

    pub fn header(&self) -> Option<&HtsHdr> {
        self.hts.header()
    }
}

