use std::{ffi::CStr, fmt, io, ops::Deref, ptr, str::FromStr};

use super::{
    check_tid, BamAux, BamAuxItem, BamAuxItemMut, BamAuxIter, BamAuxIterMut, BamRec, SamHeader,
    SeqQual,
};

use crate::{
    from_cstr, get_cstr, htsFile, hts_err, hts_idx_t, hts_itr_t, kstring_t, HtsFile, HtsPos,
};

use libc::{c_char, c_int, c_void, size_t};

pub const BAM_FPAIRED: u16 = 1;
pub const BAM_FPROPER_PAIR: u16 = 2;
pub const BAM_FUNMAP: u16 = 4;
pub const BAM_FMUNMAP: u16 = 8;
pub const BAM_FREVERSE: u16 = 16;
pub const BAM_FMREVERSE: u16 = 32;
pub const BAM_FREAD1: u16 = 64;
pub const BAM_FREAD2: u16 = 128;
pub const BAM_FSECONDARY: u16 = 256;
pub const BAM_FQCFAIL: u16 = 512;
pub const BAM_FDUP: u16 = 1024;
pub const BAM_FSUPPLEMENTARY: u16 = 2048;

pub const SAM_FORMAT_VERSION: &str = "1.6";

#[link(name = "hts")]
extern "C" {
    pub(super) fn sam_hdr_read(fp_: *mut htsFile) -> *mut sam_hdr_t;
    pub(super) fn sam_hdr_write(fp_: *mut htsFile, hd_: *const sam_hdr_t) -> c_int;
    pub(super) fn sam_hdr_init() -> *mut sam_hdr_t;
    pub(super) fn sam_hdr_destroy(hd_: *mut sam_hdr_t);
    pub(super) fn sam_hdr_dup(hd_: *const sam_hdr_t) -> *mut sam_hdr_t;
    pub(super) fn sam_hdr_add_lines(
        hd_: *mut sam_hdr_t,
        lines_: *const c_char,
        len_: size_t,
    ) -> c_int;
    pub(super) fn sam_hdr_remove_except(
        hd_: *mut sam_hdr_t,
        type_: *const c_char,
        id_key_: *const c_char,
        id_value_: *const c_char,
    ) -> c_int;
    pub(super) fn sam_hdr_nref(hd_: *const sam_hdr_t) -> c_int;
    pub(super) fn sam_hdr_tid2name(hd_: *const sam_hdr_t, i_: c_int) -> *const c_char;
    pub(super) fn sam_hdr_tid2len(hd_: *const sam_hdr_t, i_: c_int) -> c_int;
    pub(super) fn sam_hdr_name2tid(hd_: *const sam_hdr_t, nm_: *const c_char) -> c_int;
    pub(super) fn sam_hdr_str(hd_: *const sam_hdr_t) -> *const c_char;
    pub(super) fn bam_init1() -> *mut bam1_t;
    pub(super) fn bam_destroy1(b: *mut bam1_t);
    pub(super) fn bam_copy1(bdest: *mut bam1_t, bsrc: *const bam1_t) -> *mut bam1_t;
    pub(super) fn bam_endpos(pt_: *const bam1_t) -> HtsPos;
    pub(super) fn bam_aux_update_str(
        pt_: *mut bam1_t,
        tag_: *const c_char,
        len_: c_int,
        data_: *const c_char,
    ) -> c_int;
    pub(super) fn bam_aux_append(
        pt: *mut bam1_t,
        tag: *const c_char,
        type_: c_char,
        len: c_int,
        data: *const u8,
    ) -> c_int;
    pub(super) fn sam_read1(fp_: *mut htsFile, hd_: *mut sam_hdr_t, b_: *mut bam1_t) -> c_int;
    pub(super) fn sam_write1(fp_: *mut htsFile, hd_: *mut sam_hdr_t, b_: *const bam1_t) -> c_int;
    pub(super) fn sam_format1(hdr: *const sam_hdr_t, b: *const bam1_t, s: *mut kstring_t) -> c_int;
    pub(super) fn sam_open_mode(mode: *mut c_char, fn_: *const c_char, fmt: *const c_char)
        -> c_int;
    pub(super) fn sam_hdr_change_HD(hd: *mut sam_hdr_t, key: *const c_char, val: *const c_char);
    pub(super) fn sam_hdr_find_line_id(
        hd: *mut sam_hdr_t,
        type_: *const c_char,
        id_key: *const c_char,
        id_val: *const c_char,
        ks: *mut kstring_t,
    ) -> c_int;
    pub(super) fn sam_hdr_find_line_pos(
        hd: *mut sam_hdr_t,
        type_: *const c_char,
        pos: c_int,
        ks: *mut kstring_t,
    ) -> c_int;
    pub(super) fn sam_hdr_find_tag_id(
        hd: *mut sam_hdr_t,
        type_: *const c_char,
        id_key: *const c_char,
        id_val: *const c_char,
        key: *const c_char,
        ks: *mut kstring_t,
    ) -> c_int;
    pub(super) fn sam_hdr_find_tag_pos(
        hd: *mut sam_hdr_t,
        type_: *const c_char,
        pos: c_int,
        key: *const c_char,
        ks: *mut kstring_t,
    ) -> c_int;
    pub(super) fn sam_hdr_count_lines(hd: *mut sam_hdr_t, type_: *const c_char) -> c_int;
    pub(crate) fn sam_itr_queryi(
        idx: *const hts_idx_t,
        tid: c_int,
        start: HtsPos,
        end: HtsPos,
    ) -> *mut hts_itr_t;
    pub(crate) fn sam_idx_init(
        fp: *mut htsFile,
        hd: *const sam_hdr_t,
        min_shift: c_int,
        fnidx: *const c_char,
    ) -> c_int;
    pub(crate) fn sam_idx_save(fp: *mut htsFile) -> c_int;
    pub fn sam_hdr_add_pg(hd: *mut sam_hdr_t, name: *const c_char, ...) -> c_int;
}

pub(crate) unsafe extern "C" fn bam_name2id(hdr: *mut c_void, s: *const c_char) -> c_int {
    sam_hdr_name2tid(hdr as *mut sam_hdr_t, s)
}

#[repr(C)]
pub struct sam_hdr_t {
    _unused: [u8; 0],
}

impl sam_hdr_t {
    pub fn write(&self, hts_file: &mut HtsFile) -> io::Result<()> {
        match unsafe { sam_hdr_write(hts_file.as_mut(), self) } {
            0 => Ok(()),
            _ => Err(hts_err(format!(
                "Failed to write header to {}",
                hts_file.name()
            ))),
        }
    }
    pub fn nref(&self) -> usize {
        let l = unsafe { sam_hdr_nref(self) };
        l as usize
    }
    fn check_idx(&self, i: usize) {
        if i >= self.nref() {
            panic!("Reference ID {} out of range", i);
        }
    }
    pub fn tid2name(&self, i: usize) -> &str {
        self.check_idx(i);
        from_cstr(unsafe { sam_hdr_tid2name(self, i as c_int) })
    }
    pub fn tid2len(&self, i: usize) -> usize {
        self.check_idx(i);
        let len = unsafe { sam_hdr_tid2len(self, i as c_int) };
        len as usize
    }
    pub fn name2tid<S: AsRef<str>>(&self, cname: S) -> Option<usize> {
        let tid = unsafe { sam_hdr_name2tid(self, get_cstr(cname).as_ptr()) };
        if tid < 0 {
            None
        } else {
            Some(tid as usize)
        }
    }
    pub fn text(&self) -> &str {
        from_cstr(unsafe { sam_hdr_str(self) })
    }

    pub fn add_lines<S: AsRef<str>>(&mut self, lines: S) -> io::Result<()> {
        let lines = lines.as_ref();
        match unsafe { sam_hdr_add_lines(self, get_cstr(lines).as_ptr(), lines.len() as size_t) } {
            0 => Ok(()),
            _ => Err(hts_err("Failed to add line to SAM/BAM header".to_string())),
        }
    }
    pub fn remove_except(
        &mut self,
        ln_type: &str,
        id_key: Option<&str>,
        id_value: Option<&str>,
    ) -> io::Result<()> {
        match if let (Some(key), Some(value)) = (id_key, id_value) {
            unsafe {
                sam_hdr_remove_except(
                    self,
                    get_cstr(ln_type).as_ptr(),
                    get_cstr(key).as_ptr(),
                    get_cstr(value).as_ptr(),
                )
            }
        } else {
            unsafe {
                sam_hdr_remove_except(self, get_cstr(ln_type).as_ptr(), ptr::null(), ptr::null())
            }
        } {
            0 => Ok(()),
            _ => Err(hts_err(format!(
                "Failed to remove {} lines from SAM/BAM header",
                ln_type
            ))),
        }
    }
    pub fn remove(&mut self, ln_type: &str) -> io::Result<()> {
        self.remove_except(ln_type, None, None)
    }
    pub fn change_hd(&mut self, key: &str, val: Option<&str>) {
        let val = if let Some(v) = val {
            get_cstr(v).as_ptr()
        } else {
            std::ptr::null::<c_char>()
        };
        unsafe { sam_hdr_change_HD(self, get_cstr(key).as_ptr(), val) }
    }
    pub fn find_line_id(&mut self, typ: &str, id_key: &str, id_val: &str) -> Option<kstring_t> {
        let mut ks = kstring_t::new();
        if unsafe {
            sam_hdr_find_line_id(
                self,
                get_cstr(typ).as_ptr(),
                get_cstr(id_key).as_ptr(),
                get_cstr(id_val).as_ptr(),
                &mut ks,
            ) == 0
        } {
            Some(ks)
        } else {
            None
        }
    }
    pub fn find_line_pos(&mut self, typ: &str, pos: usize) -> Option<kstring_t> {
        let mut ks = kstring_t::new();
        if unsafe {
            sam_hdr_find_line_pos(self, get_cstr(typ).as_ptr(), pos as c_int, &mut ks) == 0
        } {
            Some(ks)
        } else {
            None
        }
    }
    pub fn find_tag_id(
        &mut self,
        typ: &str,
        id_key: &str,
        id_val: &str,
        key: &str,
    ) -> Option<kstring_t> {
        let mut ks = kstring_t::new();
        if unsafe {
            sam_hdr_find_tag_id(
                self,
                get_cstr(typ).as_ptr(),
                get_cstr(id_key).as_ptr(),
                get_cstr(id_val).as_ptr(),
                get_cstr(key).as_ptr(),
                &mut ks,
            ) == 0
        } {
            Some(ks)
        } else {
            None
        }
    }
    pub fn find_tag_pos(&mut self, typ: &str, pos: usize, key: &str) -> Option<kstring_t> {
        let mut ks = kstring_t::new();
        if unsafe {
            sam_hdr_find_tag_pos(
                self,
                get_cstr(typ).as_ptr(),
                pos as c_int,
                get_cstr(key).as_ptr(),
                &mut ks,
            ) == 0
        } {
            Some(ks)
        } else {
            None
        }
    }
    pub fn count_lines(&mut self, typ: &str) -> Option<usize> {
        let n = unsafe { sam_hdr_count_lines(self, get_cstr(typ).as_ptr()) };
        if n >= 0 {
            Some(n as usize)
        } else {
            None
        }
    }
}

#[repr(C)]
pub struct bam1_core_t {
    pos: HtsPos,
    tid: i32,
    bin: u16,
    qual: u8,
    l_extranul: u8,
    flag: u16,
    l_qname: u16,
    n_cigar: u32,
    l_qseq: i32,
    mtid: i32,
    mpos: HtsPos,
    isze: HtsPos,
}

#[repr(C)]
pub struct bam1_t {
    core: bam1_core_t,
    id: u64,
    pub(super) data: *mut c_char,
    l_data: c_int,
    m_data: u32,
    mempolicy: u32,
}

/// Rust representation of all fields in bam1_t data as mutable slices
/// Note: the fields are public so that they can be accessed without function calls.
/// In this way it is possible to use all 4 mutable references at the same time (this is
/// safe as we assure that they do not overlap).
pub struct BamData<'a> {
    pub qname: &'a mut [c_char],
    pub cigar: &'a mut [CigarElem],
    pub seq: &'a mut [u8],
    pub qual: &'a mut [u8],
    pub aux: &'a mut [u8],
}

impl bam1_t {
    /// Get all fields in bam1_t data as rust mutable slices allowing them to be
    /// (more or less) safely manipulated.  This can't be done using the individual
    /// function calls (i.e., bam1_t::get_seq_mut()) as the borrow checker will not allow more
    /// than one mutable reference to be held from the bam1_t structure at a time.  Here we can
    /// assure that the 5 mutable slices are non-overlapping so there is no risk of aliasing
    pub fn bam_data(&mut self) -> Option<BamData> {
        let mut p = self.data;
        if p.is_null() {
            None
        } else {
            let size = self.l_data as usize;
            unsafe {
                let end_p = p.add(size);
                let qlen = libc::strlen(p) as usize;
                let qname = std::slice::from_raw_parts_mut(p, qlen + 1);
                p = p.add(self.core.l_qname as usize);
                assert_eq!(
                    p.align_offset(4),
                    0,
                    "Cigar offset not aligned - Bam record corrupt"
                );
                let cigar =
                    std::slice::from_raw_parts_mut(p as *mut CigarElem, self.core.n_cigar as usize);
                p = p.add((self.core.n_cigar << 2) as usize);
                let qlen = self.core.l_qseq as usize;
                let sqlen = (qlen + 1) >> 1;
                let seq = std::slice::from_raw_parts_mut(p as *mut u8, sqlen);
                p = p.add(sqlen);
                let qual = std::slice::from_raw_parts_mut(p as *mut u8, qlen);
                p = p.add(qlen);
                let aux_len = end_p.offset_from(p);
                assert!(aux_len >= 0, "Corrupt BAM record");
                let aux = std::slice::from_raw_parts_mut(p as *mut u8, aux_len as usize);
                Some(BamData {
                    qname,
                    cigar,
                    seq,
                    qual,
                    aux,
                })
            }
        }
    }

    pub fn qname_cstr(&self) -> Option<&CStr> {
        if self.data.is_null() {
            None
        } else {
            Some(unsafe { CStr::from_ptr(self.data) })
        }
    }

    pub fn qname(&self) -> io::Result<&str> {
        if self.data.is_null() {
            Err(hts_err("Empty BamRec".to_string()))
        } else {
            Ok(from_cstr(self.data))
        }
    }

    pub fn endpos(&self) -> usize {
        unsafe { bam_endpos(self) as usize }
    }

    pub fn tid(&self) -> Option<usize> {
        check_tid(self.core.tid)
    }

    pub fn set_tid(&mut self, tid: usize) {
        self.core.tid = tid as c_int
    }

    pub fn mtid(&self) -> Option<usize> {
        check_tid(self.core.mtid)
    }

    pub fn qual(&self) -> u8 {
        self.core.qual
    }

    pub fn flag(&self) -> u16 {
        self.core.flag
    }

    pub fn pos(&self) -> Option<usize> {
        if self.core.pos >= 0 {
            Some(self.core.pos as usize)
        } else {
            None
        }
    }

    pub fn set_pos(&mut self, pos: usize) {
        self.core.pos = pos as HtsPos
    }

    pub fn mpos(&self) -> Option<usize> {
        if self.core.mpos >= 0 {
            Some(self.core.mpos as usize)
        } else {
            None
        }
    }

    pub fn set_mpos(&mut self, mpos: usize) {
        self.core.mpos = mpos as HtsPos
    }

    pub fn template_len(&self) -> isize {
        self.core.isze as isize
    }

    pub fn l_qseq(&self) -> i32 {
        self.core.l_qseq
    }

    pub fn qnames_eq(&self, b: &BamRec) -> io::Result<bool> {
        if self.data.is_null() || b.data.is_null() {
            Err(hts_err("Attempt to compare empty Bam Records".to_string()))
        } else {
            Ok(unsafe { libc::strcmp(self.data, b.data) } == 0)
        }
    }

    pub fn get_seq(&self) -> Option<&[u8]> {
        if self.data.is_null() {
            None
        } else {
            unsafe {
                let core = &self.core;
                let off = ((core.n_cigar as isize) << 2) + (core.l_qname as isize) as isize;
                let p = self.data.offset(off) as *const u8;
                let size = (core.l_qseq + 1) >> 1;
                Some(std::slice::from_raw_parts(p, size as usize))
            }
        }
    }

    pub fn get_seq_mut(&mut self) -> Option<&mut [u8]> {
        if self.data.is_null() {
            None
        } else {
            unsafe {
                let core = &self.core;
                let off = ((core.n_cigar as isize) << 2) + (core.l_qname as isize) as isize;
                let p = self.data.offset(off) as *mut u8;
                let size = (core.l_qseq + 1) >> 1;
                Some(std::slice::from_raw_parts_mut(p, size as usize))
            }
        }
    }

    pub fn get_qual(&self) -> Option<&[u8]> {
        if self.data.is_null() {
            None
        } else {
            unsafe {
                let core = &self.core;
                let off = ((core.n_cigar as isize) << 2)
                    + (core.l_qname as isize)
                    + ((core.l_qseq + 1) >> 1) as isize;
                let p = self.data.offset(off) as *const u8;
                let size = core.l_qseq;
                Some(std::slice::from_raw_parts(p, size as usize))
            }
        }
    }

    pub fn get_qual_mut(&mut self) -> Option<&mut [u8]> {
        if self.data.is_null() {
            None
        } else {
            unsafe {
                let core = &self.core;
                let off = ((core.n_cigar as isize) << 2)
                    + (core.l_qname as isize)
                    + ((core.l_qseq + 1) >> 1) as isize;
                let p = self.data.offset(off) as *mut u8;
                let size = core.l_qseq;
                Some(std::slice::from_raw_parts_mut(p, size as usize))
            }
        }
    }

    pub fn cigar(&self) -> Option<Cigar> {
        let len = self.core.n_cigar as usize;
        if len > 0 {
            assert!(!self.data.is_null());
            let slice = unsafe {
                let ptr = self.data.offset(self.core.l_qname as isize);
                assert_eq!(
                    ptr.align_offset(4),
                    0,
                    "Cigar storage not aligned - Bam record corrupt"
                );
                std::slice::from_raw_parts(ptr.cast::<CigarElem>(), len)
            };
            Some(Cigar(slice))
        } else {
            None
        }
    }

    pub fn qlen(&self) -> Option<u32> {
        self.cigar().map(|c| c.qlen())
    }

    pub fn rlen(&self) -> Option<u32> {
        self.cigar().map(|c| c.rlen())
    }

    pub fn cigar_buf(&self) -> Option<CigarBuf> {
        self.cigar().map(|c| c.to_cigar_buf())
    }

    pub fn write(&mut self, hfile: &mut HtsFile, hdr: &mut SamHeader) -> io::Result<usize> {
        match unsafe { sam_write1(hfile.as_mut(), hdr.as_mut(), self) } {
            x if x >= 0 => Ok(x as usize),
            _ => Err(hts_err(format!(
                "Failed to write BamRec to {}",
                hfile.name()
            ))),
        }
    }

    pub fn aux_update_str<S: AsRef<str>>(&mut self, tag: &str, data: S) -> io::Result<()> {
        if tag.len() != 2 {
            return Err(hts_err(
                "Failed to update string tag: tag length is not 2".to_string(),
            ));
        }
        let data = data.as_ref();
        match unsafe {
            bam_aux_update_str(
                self,
                get_cstr(tag).as_ptr(),
                data.len() as c_int,
                get_cstr(data).as_ptr(),
            )
        } {
            0 => Ok(()),
            _ => Err(hts_err("Failed to update string tag".to_string())),
        }
    }

    pub fn get_aux_as_parts(&self) -> Option<(*mut u8, usize)> {
        if self.data.is_null() {
            None
        } else {
            unsafe {
                let core = &self.core;
                let off = ((core.n_cigar as isize) << 2)
                    + (core.l_qname as isize)
                    + (core.l_qseq + ((core.l_qseq + 1) >> 1)) as isize;
                let p = self.data.offset(off) as *mut u8;
                let size = self.l_data as isize - off;
                assert!(size >= 0, "Invalid BAM aux size");
                Some((p, size as usize))
            }
        }
    }

    pub fn get_aux(&self) -> Option<&[u8]> {
        self.get_aux_as_parts()
            .map(|(p, sz)| unsafe { std::slice::from_raw_parts(p, sz) })
    }

    pub fn get_aux_iter(&self) -> Option<BamAuxIter> {
        self.get_aux().map(|aux| BamAuxIter { data: aux })
    }

    pub fn get_aux_iter_mut(&mut self) -> Option<BamAuxIterMut> {
        self.get_aux_as_parts()
            .and_then(|(p, sz)| BamAuxIterMut::new(p, sz))
    }

    pub fn get_tag(&self, tag_id: &str) -> Option<BamAuxItem> {
        self.get_aux_iter().and_then(|a| find_tag(a, tag_id))
    }

    pub fn get_tag_mut(&mut self, tag_id: &str) -> Option<BamAuxItemMut> {
        self.get_aux_iter_mut().and_then(|a| find_tag(a, tag_id))
    }

    /// Delete Aux tag from bam record where we already have the pointer to the start of the tag
    /// entry within the record and the size of the entry.  It is much safer to use del_tag(), but
    /// this function avoids searching through the tags if the position is already known.
    ///
    /// # Safety
    ///
    /// The supplied ptr and sz arguments must correspond to the start and size of an aux tag in self
    /// otherwise behaviour is undefined (although there are many sanity checks in place. i.e., we
    /// check that the given region lies entirely within the aux data area, but we don't check that
    /// p and sz correspond to a true starting position and size of a tag).
    pub unsafe fn del_tag_from_ptr(&mut self, p: *mut u8, sz: usize) {
        let (aux_ptr, size) = self.get_aux_as_parts().expect("Empty bam rec");
        let size = size as isize;
        let isz = sz as isize;
        let x = p.offset_from(aux_ptr) + isz;
        assert!(x <= size && x >= isz);
        if x < size {
            libc::memmove(
                p as *mut c_void,
                p.add(sz) as *const c_void,
                (size - x) as size_t,
            );
        }
        self.l_data -= sz as c_int;
    }

    pub fn del_tag(&mut self, tag_id: &str) -> io::Result<bool> {
        if let Some(a) = self.get_aux_iter_mut() {
            if let Some((p, sz)) = find_tag(a, tag_id).map(|mut t| (t.as_mut_ptr(), t.len())) {
                unsafe { self.del_tag_from_ptr(p, sz) }
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub fn add_tag(&mut self, tag_id: &str, tag_type: char, data: &[u8]) -> io::Result<()> {
        if tag_id.len() != 2 {
            Err(hts_err("Bad tag length in add_tag()".to_string()))
        } else {
            let r = unsafe {
                bam_aux_append(
                    self,
                    tag_id.as_ptr() as *const c_char,
                    tag_type as c_char,
                    data.len() as c_int,
                    data.as_ptr(),
                )
            };
            if r == 0 {
                Ok(())
            } else {
                Err(hts_err("Out of memory in add_tag()".to_string()))
            }
        }
    }

    pub fn get_seq_qual(&self) -> io::Result<SeqQual> {
        let seq = self
            .get_seq()
            .ok_or_else(|| hts_err("No Sequence data in BAM record".to_string()))?;
        let qual = self
            .get_qual()
            .ok_or_else(|| hts_err("No Quality data in BAM record".to_string()))?;
        let mut sq = Vec::with_capacity(qual.len());
        let mut qitr = qual.iter();
        for s in seq.iter() {
            let (b, a) = SEQ_DECODE[*s as usize];
            let q = (*qitr.next().unwrap()).min(62);
            sq.push(if a > 0 { (a - 1) | (q << 2) } else { 0 });
            if let Some(x) = qitr.next() {
                let q = (*x).min(62);
                sq.push(if b > 0 { (b - 1) | (q << 2) } else { 0 });
            }
        }
        Ok(SeqQual(sq.into_boxed_slice()))
    }

    pub fn format(&mut self, hdr: &SamHeader, s: &mut kstring_t) -> c_int {
        unsafe { sam_format1(hdr.as_ref(), self, s) }
    }
}

fn find_tag<'a, I, T>(a: I, s: &str) -> Option<T>
where
    I: Iterator<Item = T>,
    T: BamAux<'a>,
{
    let tag_id = s.as_bytes();
    if tag_id.len() != 2 {
        return None;
    }
    for tag in a {
        let tg = tag.tag();
        if tg[0] == tag_id[0] && tg[1] == tag_id[1] {
            return Some(tag);
        }
    }
    None
}

const SEQ_DECODE: [(u8, u8); 256] = [
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 1),
    (1, 1),
    (2, 1),
    (0, 1),
    (3, 1),
    (0, 1),
    (0, 1),
    (0, 1),
    (4, 1),
    (0, 1),
    (0, 1),
    (0, 1),
    (0, 1),
    (0, 1),
    (0, 1),
    (0, 1),
    (0, 2),
    (1, 2),
    (2, 2),
    (0, 2),
    (3, 2),
    (0, 2),
    (0, 2),
    (0, 2),
    (4, 2),
    (0, 2),
    (0, 2),
    (0, 2),
    (0, 2),
    (0, 2),
    (0, 2),
    (0, 2),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 3),
    (1, 3),
    (2, 3),
    (0, 3),
    (3, 3),
    (0, 3),
    (0, 3),
    (0, 3),
    (4, 3),
    (0, 3),
    (0, 3),
    (0, 3),
    (0, 3),
    (0, 3),
    (0, 3),
    (0, 3),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 4),
    (1, 4),
    (2, 4),
    (0, 4),
    (3, 4),
    (0, 4),
    (0, 4),
    (0, 4),
    (4, 4),
    (0, 4),
    (0, 4),
    (0, 4),
    (0, 4),
    (0, 4),
    (0, 4),
    (0, 4),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 0),
    (3, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (4, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0),
];

#[repr(u8)]
#[derive(Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match,
    Ins,
    Del,
    RefSkip,
    SoftClip,
    HardClip,
    Pad,
    Equal,
    Diff,
    Back,
    Overlap,
    Invalid1,
    Invalid2,
    Invalid3,
    Invalid4,
    Invalid5,
}

impl fmt::Display for CigarOp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                CigarOp::Match => 'M',
                CigarOp::Ins => 'I',
                CigarOp::Del => 'D',
                CigarOp::RefSkip => 'N',
                CigarOp::SoftClip => 'S',
                CigarOp::HardClip => 'H',
                CigarOp::Pad => 'P',
                CigarOp::Equal => '=',
                CigarOp::Diff => 'X',
                CigarOp::Back => 'B',
                CigarOp::Overlap => 'O',
                _ => '?',
            }
        )
    }
}

#[derive(Debug, Copy, Clone)]
pub struct CigarElem(u32);

const CIGAR_TYPE: u32 = 0x13C1A7;
const CIGAR_TYPE1: u32 = 0x13C5A7;

impl CigarElem {
    pub fn op_len(&self) -> u32 {
        self.0 >> 4
    }
    pub fn op(&self) -> CigarOp {
        unsafe { std::mem::transmute((self.0 & 15) as u8) }
    }
    pub fn op_pair(&self) -> (CigarOp, u32) {
        (self.op(), self.op_len())
    }

    // This magic comes from htslib/sam.h
    // If bit 0 is set in op_type then the op consumes the query, and
    // if bit 1 is set then the op consumes the reference
    pub fn op_type(&self) -> u32 {
        (CIGAR_TYPE >> ((self.0 & 15) << 1)) & 3
    }
    // Similar to above, but we also count Hard clips the same as Soft clips
    pub fn op_type1(&self) -> u32 {
        (CIGAR_TYPE1 >> ((self.0 & 15) << 1)) & 3
    }
}

const CIGAR_STR: &str = "MIDNSHP=XDO?????";

impl fmt::Display for CigarElem {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}{}",
            self.0 >> 4,
            CIGAR_STR.as_bytes()[(self.0 & 15) as usize] as char
        )
    }
}

pub struct Cigar<'a>(&'a [CigarElem]);

impl<'a> Deref for Cigar<'a> {
    type Target = [CigarElem];
    fn deref(&self) -> &[CigarElem] {
        self.0
    }
}

impl<'a> Cigar<'a> {
    pub fn qlen(&self) -> u32 {
        self.iter()
            .filter(|c| (c.op_type() & 1) != 0)
            .fold(0, |mut l, c| {
                l += c.op_len();
                l
            })
    }
    pub fn qlen1(&self) -> u32 {
        self.iter()
            .filter(|c| (c.op_type1() & 1) != 0)
            .fold(0, |mut l, c| {
                l += c.op_len();
                l
            })
    }
    pub fn rlen(&self) -> u32 {
        self.iter()
            .filter(|c| (c.op_type() & 2) != 0)
            .fold(0, |mut l, c| {
                l += c.op_len();
                l
            })
    }
    pub fn to_cigar_buf(&self) -> CigarBuf {
        let v = self.to_vec();
        CigarBuf(v.into_boxed_slice())
    }
}

pub struct CigarBuf(Box<[CigarElem]>);

impl Deref for CigarBuf {
    type Target = [CigarElem];
    fn deref(&self) -> &[CigarElem] {
        self.0.deref()
    }
}

impl fmt::Display for CigarBuf {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for elem in self.iter() {
            write!(f, "{}", elem)?
        }
        Ok(())
    }
}

const BAM_COVERLAP: u32 = 10;

fn trim_cigar_vec<I: Iterator<Item = CigarElem>>(it: I, x: u32) -> Vec<CigarElem> {
    let mut ct = 0;
    let mut v = Vec::new();
    for elem in it {
        if ct >= x || (elem.op_type() & 2) == 0 {
            v.push(elem)
        } else {
            let l = elem.op_len();
            if (elem.op_type() & 1) != 0 {
                if ct + l <= x {
                    v.push(CigarElem((elem.0 & 0xfffffff0) | BAM_COVERLAP));
                } else {
                    v.push(CigarElem(((x - ct) << 4) | BAM_COVERLAP));
                    v.push(CigarElem(((ct + l - x) << 4) | (elem.0 & 15)));
                }
            }
            ct += l;
        }
    }
    v
}

impl CigarBuf {
    pub fn qlen(&self) -> u32 {
        self.iter()
            .filter(|c| (c.op_type() & 1) != 0)
            .fold(0, |mut l, c| {
                l += c.op_len();
                l
            })
    }
    pub fn qlen1(&self) -> u32 {
        self.iter()
            .filter(|c| (c.op_type1() & 1) != 0)
            .fold(0, |mut l, c| {
                l += c.op_len();
                l
            })
    }
    pub fn rlen(&self) -> u32 {
        self.iter()
            .filter(|c| (c.op_type() & 2) != 0)
            .fold(0, |mut l, c| {
                l += c.op_len();
                l
            })
    }

    // Adjust cigar so that alignment starts x bases later w.r.t the reference
    pub fn trim_start(&mut self, x: u32) {
        let v = trim_cigar_vec(self.iter().copied(), x);
        self.0 = v.into_boxed_slice();
    }
    // Adjust cigar so that alignment ends x bases earlier w.r.t the reference
    pub fn trim_end(&mut self, x: u32) {
        let mut v = trim_cigar_vec(self.iter().copied().rev(), x);
        let v1: Vec<_> = v.drain(..).rev().collect();
        self.0 = v1.into_boxed_slice();
    }
}

const BAM_CIGAR_TAB: [i8; 256] = [
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 7, -1, -1, -1, -1, 9, -1, 2, -1, -1, -1, 5,
    1, -1, -1, -1, 0, 3, -1, 6, -1, -1, 4, -1, -1, -1, -1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
];

const BAM_CIGAR_MAX_LEN: u32 = 1 << 28;

impl FromStr for CigarBuf {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut sp = s;
        let mut v = Vec::new();
        while !sp.is_empty() {
            let i = sp
                .find(|c: char| !c.is_ascii_digit())
                .ok_or("Cigar string does not end in letter")?;
            let n = <u32>::from_str(&sp[0..i])
                .map_err(|_| "Error parsing Cigar string - expecting number")?;
            if n >= BAM_CIGAR_MAX_LEN {
                return Err("Cigar number too large");
            };
            let op = BAM_CIGAR_TAB[sp[i..=i].as_bytes()[0] as usize];
            if op < 0 {
                return Err("Illegal Cigar character");
            }
            v.push(CigarElem((n << 4) | (op as u32)));
            sp = &sp[i + 1..];
        }
        Ok(CigarBuf(v.into_boxed_slice()))
    }
}
