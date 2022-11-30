pub mod lib;
pub use lib::*;

use libc::{c_void, c_int, size_t};

use std::{
    io:: {self, Write},
    marker::PhantomData,
    ptr::{null, NonNull},
    ops::{Deref, DerefMut},
    ffi::{CStr, CString},
    cmp::Ordering,
    path::{Path, PathBuf},
    os::unix::ffi::OsStrExt,
};

use crate::{
    tbx_index_load3, VcfHeader, SamHeader, Tbx, vcf_read_itr, bcf_read_itr, tbx_read_itr, tbx_name2id,
    sam_itr_queryi, sam_hdr_t, bam_name2id, BamRec, BcfRec, TbxRec, bcf_hdr_name2id,
    get_cstr, hts_err, bcf_hdr_t, tbx_t, BCF_DT_CTG, bcf_idx_init, bcf_idx_save, sam_idx_init, sam_idx_save
};

pub struct Hts {
    pub(crate) hts_file: HtsFile,
    pub(crate) header: Option<HtsHdr>,
    pub(crate) tbx: Option<Tbx>,
    pub(crate) idx: Option<HtsIndex>,
}

unsafe impl Send for Hts {}
unsafe impl Sync for Hts {}

impl Hts {
    pub fn open<P>(name: Option<P>, mode: &str) -> io::Result<Self>
        where
            P: AsRef<Path>,
    {
        Self::open_format_(name, mode, None)
    }

    pub fn open_format<P>(name: Option<P>, mode: &str, fmt: &HtsFormat) -> io::Result<Self>
        where
            P: AsRef<Path>,
    {
        Self::open_format_(name, mode, Some(fmt))
    }

    fn open_format_<P>(name: Option<P>, mode: &str, fmt: Option<&HtsFormat>) -> io::Result<Self>
        where
            P: AsRef<Path>,
    {
        let name = name.map(|p| p.as_ref().to_owned()).unwrap_or_else(|| PathBuf::from("-"));
        let mut fp = HtsFile::open_format_(&name, mode, fmt)?;
        let (header, tbx) = if fp.as_ref().is_write() == 0 {
            match fp.format().format {
                htsExactFormat::Sam | htsExactFormat::Bam | htsExactFormat::Cram => (Some(HtsHdr::Sam(SamHeader::read(&mut fp)?)), None),
                htsExactFormat::Bcf | htsExactFormat::Vcf => (Some(HtsHdr::Vcf(VcfHeader::read(&mut fp)?)), None),
                _ => if let Ok(tbx) = Tbx::load(&name) { (None, Some(tbx)) } else { (None, None) }
            }
        } else { (None, None) };

        Ok(Self {
            hts_file: fp,
            header,
            tbx,
            idx: None,
        })
    }

    pub fn set_header(&mut self, hdr: Option<HtsHdr>) {
        self.header = hdr
    }

    pub fn header(&self) -> Option<&HtsHdr> { self.header.as_ref() }

    pub fn header_mut(&mut self) -> Option<&mut HtsHdr> { self.header.as_mut() }

    pub fn seq_names(&self) -> Vec<&str> {
        if let Some(h) = self.header.as_ref() { h.seq_names() }
        else if let Some(t) = self.tbx.as_ref() { t.seq_names() }
        else { Vec::new() }
    }

    pub fn seq_lengths(&self) -> Vec<usize> {
        if let Some(h) = self.header.as_ref() { h.seq_lengths() }
        else { Vec::new() }
    }

    pub fn seq_length(&self, ctg: &str) -> Option<usize> {
        self.header.as_ref().and_then(|h| h.seq_length(ctg))
    }

    pub fn name2tid(&self, ctg: &str) -> Option<usize> {
        self.header.as_ref().and_then(|h| h.name2tid(ctg))
    }

    pub fn n_ref(&self) -> Option<usize> {
        self.header.as_ref().map(|h| h.n_ref())
    }

    pub fn tbx(&self) -> Option<&Tbx> { self.tbx.as_ref() }

    pub fn tbx_mut(&mut self) -> Option<&mut Tbx> { self.tbx.as_mut() }

    pub fn hts_file_mut(&mut self) -> &mut HtsFile { &mut self.hts_file }

    pub fn hts_file(&self) -> &HtsFile { &self.hts_file }

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

    pub fn idx_init(&mut self, min_shift: isize, fnidx: &CStr) -> io::Result<()> {
        let (hts_file, hdr) = self.hts_file_and_header();
        if let Some(h) = hdr {
            let ret = match h {
                HtsHdr::Vcf(h) => unsafe { bcf_idx_init(hts_file.as_mut(), h.as_ref(), min_shift as c_int, fnidx.as_ptr())}
                HtsHdr::Sam(h) => unsafe { sam_idx_init(hts_file.as_mut(), h.as_ref(), min_shift as c_int, fnidx.as_ptr())}
            };
            if ret == 0 {
                return Ok(())
            }
        }
        Err(hts_err(format!("Failed to initialize index {:?} - wrong file type", fnidx)))
    }

    pub fn idx_save(&mut self) -> io::Result<()> {
        let hts_file = &mut self.hts_file;
        let ret = match hts_file.format().format {
            htsExactFormat::Vcf | htsExactFormat::Bcf => unsafe { bcf_idx_save(hts_file.as_mut()) },
            _ => unsafe {sam_idx_save(hts_file.as_mut())}
        };
        if ret == 0 {
            Ok(())
        } else {
            Err(hts_err("Failed to save index".to_string()))
        }
    }

    pub fn rec_type(&self) -> Option<HtsRecType> {
        if self.tbx.is_some() { Some(HtsRecType::Tbx) } else {
            self.header().map(|h| match h {
                HtsHdr::Sam(_) => HtsRecType::Sam,
                HtsHdr::Vcf(_) => HtsRecType::Vcf,
            })
        }
    }

    pub fn itr_query(&mut self, region: &Region) -> io::Result<HtsItr> {
        if !self.has_index() { self.index_load()? }
        let idx= self.index().unwrap();
        let hdr = self.header.as_ref();

        if let Some(itr) = NonNull::new(match (hdr, self.tbx.as_ref()) {
            (_, Some(_)) => unsafe { hts_itr_query(idx, region.tid, region.begin, region.end, tbx_read_itr) },
            (Some(HtsHdr::Sam(_)), _) => unsafe { sam_itr_queryi(idx, region.tid, region.begin, region.end) },
            (Some(HtsHdr::Vcf(_)), _) if matches!(self.hts_file.format().format, htsExactFormat::Bcf) => unsafe { hts_itr_query(idx, region.tid, region.begin, region.end, bcf_read_itr) },
            (Some(HtsHdr::Vcf(_)), _) => unsafe { hts_itr_query(idx, region.tid, region.begin, region.end, vcf_read_itr) },
            (None, None) => return Err(hts_err(format!("Error making iterator for file {}", self.name()))),
        }) {
            Ok(HtsItr { inner: itr, phantom: PhantomData })
        } else {
            Err(hts_err(format!("Error making iterator for file {}", self.name())))
        }
    }

    pub fn itr_querys(&mut self, reg: &CStr) -> io::Result<HtsItr> {
        match reg.to_bytes() {
            [b'.'] => {
                let region = Region::make(HTS_IDX_START, 0, 0);
                self.itr_query(&region)
            },
            [b'*'] => {
                let region = Region::make(HTS_IDX_NOCOOR, 0, 0);
                self.itr_query(&region)
            },
            _ => {
                let region = self.parse_region(reg)?;
                self.itr_query(&region)
            },
        }
    }

    fn get_parse_data(&mut self) -> io::Result<(HtsName2Id, *mut c_void)> {
        Ok(if let Some(tbx) = self.tbx.as_mut() {
            // Tabix file
            (tbx_name2id, tbx.as_mut() as *mut tbx_t as *mut c_void)
        } else {
            // VCF/BCF or SAM/BAM/CRAM files
            match self.header_mut() {
                Some(HtsHdr::Vcf(h)) => (bcf_hdr_name2id, h.as_mut() as *mut bcf_hdr_t as *mut c_void),
                Some(HtsHdr::Sam(h)) => (bam_name2id, h.as_mut() as *mut sam_hdr_t as *mut c_void),
                _ => return Err(hts_err(format!("Error making iterator for file {}", self.name()))),
            }
        })
    }

    fn _parse_region(get_id: HtsName2Id, hdr: *mut c_void, reg: &CStr) -> io::Result<Region> {
        let mut region = Region::new();

        if unsafe {
            hts_parse_region(reg.as_ptr(), &mut region.tid, &mut region.begin, &mut region.end, get_id, hdr, HTS_PARSE_THOUSANDS_SEP)
        }.is_null() { return Err(hts_err(format!("Error parsing region {:?}", reg))) }

        Ok(region)
    }

    fn parse_region(&mut self, reg: &CStr) -> io::Result<Region> {
        let (get_id, hdr) = self.get_parse_data()?;
        Self::_parse_region(get_id, hdr, reg)
    }

    pub fn make_region_list<S: AsRef<str>>(&mut self, regs: &[S]) -> RegionList {
        let mut rlist = RegionList::new();
        if let Ok((get_id, hdr)) = self.get_parse_data() {
            for (ix, reg) in regs.iter().enumerate() {
                let r = get_cstr(reg);
                if let Ok(mut region) = match r.to_bytes() {
                    [b'.'] => Ok(Region::make(HTS_IDX_START, 0, 0)),
                    [b'*'] => Ok(Region::make(HTS_IDX_NOCOOR, 0, 0)),
                    _ => Self::_parse_region(get_id, hdr, r.as_ref()),
                } {
                    region.set_idx(Some(ix as u32));
                    rlist.add_region(region)
                }
            }
        }
        rlist
    }

    pub fn reader<T: HtsRead>(&mut self) -> HtsReader<T> { HtsReader::new(self) }

    pub fn itr_reader<'a, 'b, T>(&'a mut self, regions: &'b [Region]) -> HtsItrReader<'a, 'b, T> { HtsItrReader::new(self, regions) }

    pub fn writer(&mut self) -> io::Result<Writer> { Writer::new(&mut self.hts_file) }
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
unsafe impl Sync for HtsFile {}

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
        unsafe {
            if self.as_ref().is_write() == 0 {
                hts_flush(self.as_mut());
            }
            hts_close(self.as_mut())
        };
    }
}

impl HtsFile {
    pub fn open<S: AsRef<Path>>(name: S, mode: &str) -> io::Result<Self> {
        Self::open_format_(name, mode, None)
    }

    pub fn open_format<S: AsRef<Path>>(name: S, mode: &str, fmt: &HtsFormat) -> io::Result<Self> {
        Self::open_format_(name, mode, Some(fmt))
    }

    fn open_format_<S: AsRef<Path>>(name: S, mode: &str, fmt: Option<&HtsFormat>) -> io::Result<Self> {
        let name = name.as_ref();
        let mode = get_cstr(mode);
        let fmt = match fmt {
            Some(f) => f.inner.as_ptr(),
            None => null(),
        };
        let cname = CString::new(name.as_os_str().as_bytes()).unwrap();
        if let Some(hfile) = NonNull::new(unsafe { hts_open_format(cname.as_ptr(), mode.as_ptr(), fmt) }) {
            assert!(!unsafe {hfile.as_ref()}.fn_.is_null());
            let hfile = hfile.cast::<htsFile>();
            Ok(HtsFile {
                inner: hfile,
                phantom: PhantomData,
            })
        } else {
            Err(hts_err(format!("Couldn't open file {}", name.display())))
        }
    }

    pub fn set_thread_pool(&mut self, tp: &HtsThreadPool) -> io::Result<()> {
        if unsafe{hts_set_thread_pool(self.as_mut(), tp)} == 0 {
            Ok(())
        } else {
            Err(hts_err(format!("Failed to set thread pool for file {}", self.name())))
        }
    }

    pub fn writer(&mut self) -> io::Result<Writer> { Writer::new(self) }
}

#[derive(Debug, Clone)]
pub enum HtsHdr {
    Vcf(VcfHeader),
    Sam(SamHeader),
}

impl HtsHdr {
    pub fn n_ref(&self) -> usize {
        match self {
            HtsHdr::Vcf(h) => h.n_ref() as usize,
            HtsHdr::Sam(h) => h.nref(),
        }
    }

    pub fn seq_names(&self) -> Vec<&str> {
        match self {
            HtsHdr::Vcf(h) => h.seq_names(),
            HtsHdr::Sam(h) => h.seq_names(),
        }
    }

    pub fn seq_lengths(&self) -> Vec<usize> {
        match self {
            HtsHdr::Vcf(h) => h.seq_lengths(),
            HtsHdr::Sam(h) => h.seq_lengths(),
        }
    }

    pub fn seq_length(&self, ctg: &str) -> Option<usize> {
        match self {
            HtsHdr::Vcf(h) => h.id2int(BCF_DT_CTG as usize, ctg).map(|i| h.id2len(i)),
            HtsHdr::Sam(h) => h.name2tid(ctg).map(|i| h.tid2len(i)),
        }
    }

    pub fn name2tid(&self, ctg: &str) -> Option<usize> {
        match self {
            HtsHdr::Vcf(h) => h.id2int(BCF_DT_CTG as usize, ctg),
            HtsHdr::Sam(h) => h.name2tid(ctg),
        }
    }
}

#[repr(C)]
#[derive(Debug)]
pub struct HtsThreadPool {
    inner: NonNull<hts_tpool>,
    qsize: c_int,
    phantom: PhantomData<hts_tpool>,
}

unsafe impl Send for HtsThreadPool {}
unsafe impl Sync for HtsThreadPool {}

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
unsafe impl Sync for HtsIndex {}

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

pub struct HtsItr {
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

pub trait HtsWrite {
    fn write(&mut self, hts: &mut Hts) -> io::Result<()>;
}

pub struct HtsReader<'a, T> {
    hts: &'a mut Hts,
    _phantom: PhantomData<T>,
}

pub trait HtsIterator<T> {
    fn read(&mut self, rec: &mut T) -> io::Result<bool>;
    fn header(&self) -> Option<&HtsHdr>;
}

impl <'a, T: HtsRead> HtsIterator<T> for HtsReader<'a, T> {
    fn read(&mut self, rec: &mut T) -> io::Result<bool> {
        rec.read(self.hts)
    }

    fn header(&self) -> Option<&HtsHdr> {
        self.hts.header()
    }
}

impl <'a, T: HtsRead> HtsReader<'a, T> {
    pub fn new(hts: &'a mut Hts) -> Self {
        Self {
            hts,
            _phantom: PhantomData
        }
    }
}

pub struct HtsItrReader<'a, 'b, T>
{
    hts: &'a mut Hts,
    region_itr: std::slice::Iter<'b, Region>,
    itr: Option<HtsItr>,
    region: Option<&'b Region>,
    _phantom: PhantomData<T>,
}

impl <'a, 'b, T> HtsItrReader<'a, 'b, T> {
    pub fn new(hts: &'a mut Hts, regions: &'b [Region]) -> Self
    {
        let region_itr = regions.into_iter();
        Self {
            hts,
            region_itr,
            itr: None,
            region: None,
            _phantom: PhantomData
        }
    }

    pub fn current_region(&self) -> Option<&Region> { self.region }
}

impl <'a, 'b, T: HtsRead> HtsIterator<T> for HtsItrReader<'a, 'b, T> {
    fn read(&mut self, rec: &mut T) -> io::Result<bool> {
        loop {
            if self.itr.is_none() {
                if let Some(reg) = self.region_itr.next() {
                    self.itr = Some(self.hts.itr_query(reg)?);
                    self.region = Some(reg)
                } else {
                    break
                }
            }
            let itr= self.itr.as_mut().unwrap();
            let r = rec.read_itr(self.hts, itr)?;
            if r { return Ok(true) }
            self.itr = None;
            self.region = None;
        }
        Ok(false)
    }


    fn header(&self) -> Option<&HtsHdr> {
        self.hts.header()
    }
}

#[derive(Debug, Default, Copy, Clone, Eq, PartialEq)]
pub struct Region {
    tid: c_int,
    begin: HtsPos,
    end: HtsPos,
    idx: Option<u32>,
}

impl Region {
    pub fn new() -> Self { Self::default() }

    pub fn make(tid: c_int, begin: HtsPos, end: HtsPos) -> Self {
        Self { tid, begin, end, idx: None }
    }

    pub fn set_idx(&mut self, idx: Option<u32>) { self.idx = idx }

    pub fn idx(&self) -> Option<u32> { self.idx }

    pub fn tid(&self) -> c_int { self.tid }

    pub fn begin(&self) -> HtsPos { self.begin }

    pub fn end(&self) -> HtsPos { self.end }
}

impl Ord for Region {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.tid.cmp(&other.tid) {
            Ordering::Equal => match self.begin.cmp(&other.begin) {
                Ordering::Equal => self.end.cmp(&other.end),
                ord => ord,
            },
            ord => ord,
        }
    }
}

impl PartialOrd for Region {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
}

#[derive(Debug, Default)]
pub struct RegionList {
    regions: Vec<Region>,
}

impl Deref for RegionList {
    type Target = [Region];
    #[inline]
    fn deref(&self) -> &[Region] { &self.regions }
}

impl RegionList {
    pub fn new() -> Self { Self::default() }

    pub fn add_region(&mut self, reg: Region) { self.regions.push(reg) }

    pub fn merge(&mut self) {
        if self.regions.len() > 1 {
            // Sort regions
            self.regions.sort();
            // Check for overlaps
            if self.regions.windows(2).any(|v| v[0].tid == v[1].tid && v[0].end >= v[1].begin) {
                // If we have overlaps, construct a new list, merging overlapping regions
                let mut nlist = Vec::new();
                let mut prev = self.regions[0];
                for r in self.regions.drain(1..) {
                    if prev.tid == r.tid && prev.end >= r.begin {
                        prev.end = prev.end.max(prev.end)
                    } else {
                        nlist.push(prev);
                        prev = r
                    }
                }
                nlist.push(prev);
                self.regions = nlist;
            }
            self.regions.iter_mut().for_each(|r| r.set_idx(None));
        }
    }
}

pub enum WriterFd {
    Bgzf(NonNull<BGZF>),
    Hfile(NonNull<hfile>),
}

impl WriterFd {
    fn from_htsfile(htsfile: &mut HtsFile) -> io::Result<Self> {
        if htsfile.is_write() == 0 {
            Err(hts_err("Can not set up Writer for a read only file".to_string()))
        } else {
            match htsfile.file_desc() {
                Some(HtsFileDesc::Hfile(fd)) => Ok(WriterFd::Hfile(fd)),
                Some(HtsFileDesc::Bgzf(fd)) => Ok(WriterFd::Bgzf(fd)),
                Some(_) => Err(hts_err("Bad file type for Writer".to_string())),
                None => Err(hts_err("Null file descriptor for Writer".to_string())),
            }
        }
    }

    pub fn set_mt(&mut self, n_threads: usize, n_blocks: usize) {
        if let Self::Bgzf(fd) = self {
            unsafe { bgzf_mt(fd.as_mut(), n_threads as c_int, n_blocks as c_int ); }
        }
    }
}

impl Write for WriterFd {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let i = unsafe {match self {
            WriterFd::Hfile(mut fd) => {
                hwrite(fd.as_mut(), buf.as_ptr(), buf.len() as size_t)
            },
            WriterFd::Bgzf(mut fd) => {
                bgzf_write(fd.as_mut(), buf.as_ptr() as *const c_void, buf.len() as size_t)
            },
        }};
        if i < 0 {
            Err(hts_err("Write error".to_string()))
        } else {
            Ok(i as usize)
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        if unsafe { match self {
            WriterFd::Hfile(mut fd) => hflush(fd.as_mut()),
            WriterFd::Bgzf(mut fd) => bgzf_flush(fd.as_mut()),
        }} == 0 { Ok(()) } else { Err(hts_err("flush returned error".to_string())) }
    }
}
pub struct Writer<'a> {
    fd: WriterFd,
    phantom: PhantomData<&'a mut HtsFile>,
}

impl <'a>Writer<'a> {
    pub fn new(htsfile: &'a mut HtsFile) -> io::Result<Self> {
        let fd = WriterFd::from_htsfile(htsfile)?;
        Ok(Self{fd, phantom: PhantomData})
    }

    pub fn set_mt(&mut self, n_threads: usize, n_blocks: usize) { self.fd.set_mt(n_threads,n_blocks) }
}

impl <'a> Write for Writer<'a> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> { self.fd.write(buf) }
    fn flush(&mut self) -> io::Result<()> { self.fd.flush() }
}

pub struct OwnedWriter {
    hts: Hts,
    fd: WriterFd,
}

impl OwnedWriter {
    pub fn new(mut hts: Hts) -> io::Result<Self> {
        let htsfile = hts.hts_file_mut();
        let fd = WriterFd::from_htsfile(htsfile)?;
        Ok(Self{hts, fd})
    }

    pub fn hts(&self) -> &Hts { &self.hts }

    pub fn hts_mut(&mut self) -> &mut Hts { &mut self.hts }

    pub fn set_mt(&mut self, n_threads: usize, n_blocks: usize) { self.fd.set_mt(n_threads,n_blocks) }

}

impl Write for OwnedWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> { self.fd.write(buf) }
    fn flush(&mut self) -> io::Result<()> { self.fd.flush() }
}

pub fn hts_set_log_level(level: htsLogLevel) {
    set_log_level(level)
}

pub fn hts_get_log_level() -> htsLogLevel {
    get_log_level()
}
