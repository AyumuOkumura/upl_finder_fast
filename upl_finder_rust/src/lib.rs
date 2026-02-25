use aho_corasick::{AhoCorasick, AhoCorasickBuilder};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyList, PyModule};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::fs;

#[derive(Clone, Debug)]
pub struct Probe {
    pub probe_id: String,
    pub probe_seq: String,
}

#[derive(Clone, Debug)]
struct PatternMeta {
    probe_index: usize,
    strand: char,
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct UplFinder {
    ac: Option<AhoCorasick>,
    probes: Vec<Probe>,
    pattern_meta: Vec<PatternMeta>,
}

#[pyclass]
#[derive(Clone, Debug, Serialize)]
pub struct ProbeMatch {
    pub probe_id: String,
    pub probe_seq: String,
    pub start: usize,
    pub end: usize,
    pub strand: String,
}

#[derive(Clone, Debug)]
struct Candidate {
    left_seq: String,
    right_seq: String,
    product_size: i64,
    left_len: i64,
    right_len: i64,
    tm_left: Option<f64>,
    tm_right: Option<f64>,
    gc_left: Option<f64>,
    gc_right: Option<f64>,
    pair_penalty: Option<f64>,
    left_start: i64,
}

#[derive(Clone, Debug)]
struct RankedResult {
    left_seq: String,
    right_seq: String,
    product_size: i64,
    tm_left: Option<f64>,
    tm_right: Option<f64>,
    gc_left: Option<f64>,
    gc_right: Option<f64>,
    pair_penalty: Option<f64>,
    upl_probe_id: String,
    upl_probe_seq: String,
    upl_probe_strand: String,
    probe_start_in_amplicon: i64,
    probe_end_in_amplicon: i64,
    dist_left_3p_to_probe: i64,
    dist_probe_to_right_3p: i64,
    left_start: i64,
    left_end: i64,
    right_start: i64,
    right_end: i64,
    left_len: i64,
    right_len: i64,
    amplicon_start: i64,
    amplicon_end: i64,
    junction_spanning: bool,
    junction_spanning_detail: String,
    score: f64,
}

#[pymethods]
impl ProbeMatch {
    #[getter]
    fn probe_id(&self) -> &str {
        &self.probe_id
    }

    #[getter]
    fn probe_seq(&self) -> &str {
        &self.probe_seq
    }

    #[getter]
    fn start(&self) -> usize {
        self.start
    }

    #[getter]
    fn end(&self) -> usize {
        self.end
    }

    #[getter]
    fn strand(&self) -> &str {
        &self.strand
    }
}

fn normalize_seq(seq: &str) -> String {
    let mut s = seq.to_uppercase();
    s.retain(|c| c != ' ' && c != '\n');
    s
}

fn revcomp(seq: &str) -> String {
    let mut out = String::with_capacity(seq.len());
    for base in seq.chars().rev() {
        let rc = match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => 'N',
        };
        out.push(rc);
    }
    out
}

fn build_patterns(probes: &[Probe]) -> (Vec<String>, Vec<PatternMeta>) {
    let mut patterns: Vec<String> = Vec::with_capacity(probes.len() * 2);
    let mut pattern_meta: Vec<PatternMeta> = Vec::with_capacity(probes.len() * 2);

    for (i, probe) in probes.iter().enumerate() {
        patterns.push(probe.probe_seq.clone());
        pattern_meta.push(PatternMeta {
            probe_index: i,
            strand: '+',
        });

        let rc = revcomp(&probe.probe_seq);
        if rc != probe.probe_seq {
            patterns.push(rc);
            pattern_meta.push(PatternMeta {
                probe_index: i,
                strand: '-',
            });
        }
    }
    (patterns, pattern_meta)
}

fn build_matcher(probes: &[Probe]) -> Result<(AhoCorasick, Vec<PatternMeta>), String> {
    let (patterns, pattern_meta) = build_patterns(probes);
    if patterns.is_empty() {
        return Err("no patterns".to_string());
    }
    let pattern_refs: Vec<&str> = patterns.iter().map(|s| s.as_str()).collect();
    let ac = AhoCorasickBuilder::new()
        .build(&pattern_refs)
        .map_err(|e| e.to_string())?;
    Ok((ac, pattern_meta))
}

fn collect_matches_from_normalized(
    ac: &AhoCorasick,
    probes: &[Probe],
    pattern_meta: &[PatternMeta],
    seq: &str,
) -> Vec<ProbeMatch> {
    let mut hits: Vec<ProbeMatch> = Vec::new();
    for mat in ac.find_overlapping_iter(seq) {
        let idx = mat.pattern().as_usize();
        let Some(meta) = pattern_meta.get(idx) else {
            continue;
        };
        let Some(probe) = probes.get(meta.probe_index) else {
            continue;
        };
        hits.push(ProbeMatch {
            probe_id: probe.probe_id.clone(),
            probe_seq: probe.probe_seq.clone(),
            start: mat.start(),
            end: mat.end().saturating_sub(1),
            strand: meta.strand.to_string(),
        });
    }
    hits.sort_by(|a, b| {
        a.start
            .cmp(&b.start)
            .then_with(|| a.probe_id.cmp(&b.probe_id))
            .then_with(|| a.strand.cmp(&b.strand))
    });
    hits
}

pub fn load_upl_probes(path: &str) -> Result<Vec<Probe>, Box<dyn std::error::Error>> {
    let p = std::path::Path::new(path);
    let ext = p
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if ext.is_empty() {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "probe file has no extension").into());
    }

    let probes = if ext == "json" {
        let text = fs::read_to_string(path)?;
        parse_probes_from_json_text(&text)?
    } else if ext == "tsv" || ext == "txt" {
        let text = fs::read_to_string(path)?;
        parse_probes_from_tsv_text(&text)?
    } else {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "unsupported probe file extension (supported: .json, .tsv, .txt)",
        )
        .into());
    };
    if probes.is_empty() {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "probe file contained no probes").into());
    }
    Ok(probes)
}

fn parse_probes_from_json_text(text: &str) -> Result<Vec<Probe>, Box<dyn std::error::Error>> {
    let map: HashMap<String, String> = serde_json::from_str(text)?;
    Ok(map
        .into_iter()
        .map(|(k, v)| Probe {
            probe_id: k,
            probe_seq: normalize_seq(&v),
        })
        .collect())
}

fn parse_probes_from_tsv_text(text: &str) -> Result<Vec<Probe>, Box<dyn std::error::Error>> {
    let mut out: Vec<Probe> = Vec::new();
    for (line_idx, raw_line) in text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (id, seq) = if line.contains('\t') {
            let mut parts = line.split('\t');
            (parts.next().unwrap_or("").trim(), parts.next().unwrap_or("").trim())
        } else {
            let mut parts = line.split_whitespace();
            (parts.next().unwrap_or("").trim(), parts.next().unwrap_or("").trim())
        };
        if id.is_empty() || seq.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("invalid TSV at line {}: expected <probe_id>\\t<sequence>", line_idx + 1),
            )
            .into());
        }
        out.push(Probe {
            probe_id: id.to_string(),
            probe_seq: normalize_seq(seq),
        });
    }
    Ok(out)
}

pub fn find_upl_matches(sequence: &str, probes: &[Probe]) -> Vec<ProbeMatch> {
    if probes.is_empty() {
        return Vec::new();
    }
    let Ok((ac, pattern_meta)) = build_matcher(probes) else {
        return Vec::new();
    };
    let seq = normalize_seq(sequence);
    collect_matches_from_normalized(&ac, probes, &pattern_meta, &seq)
}

fn parse_probes_from_py(probes: &Bound<'_, PyAny>) -> PyResult<Vec<Probe>> {
    let dict = probes
        .downcast::<PyDict>()
        .map_err(|_| PyValueError::new_err("probes must be a dict[str, str]"))?;
    let mut out: Vec<Probe> = Vec::with_capacity(dict.len());
    for (k, v) in dict.iter() {
        let key: String = k.extract()?;
        let val: String = v.extract()?;
        out.push(Probe {
            probe_id: key,
            probe_seq: normalize_seq(&val),
        });
    }
    Ok(out)
}

fn get_required_string(dict: &Bound<'_, PyDict>, key: &str) -> PyResult<String> {
    let value = dict
        .get_item(key)?
        .ok_or_else(|| PyValueError::new_err(format!("candidate missing key: {key}")))?;
    value.extract::<String>()
}

fn get_required_i64(dict: &Bound<'_, PyDict>, key: &str) -> PyResult<i64> {
    let value = dict
        .get_item(key)?
        .ok_or_else(|| PyValueError::new_err(format!("candidate missing key: {key}")))?;
    value.extract::<i64>()
}

fn get_optional_f64(dict: &Bound<'_, PyDict>, key: &str) -> PyResult<Option<f64>> {
    let Some(value) = dict.get_item(key)? else {
        return Ok(None);
    };
    if value.is_none() {
        return Ok(None);
    }
    Ok(Some(value.extract::<f64>()?))
}

fn parse_candidates_from_py(candidates: &Bound<'_, PyAny>) -> PyResult<Vec<Candidate>> {
    let list = candidates
        .downcast::<PyList>()
        .map_err(|_| PyValueError::new_err("candidates must be a list[dict]"))?;

    let mut out = Vec::with_capacity(list.len());
    for item in list.iter() {
        let dict = item
            .downcast::<PyDict>()
            .map_err(|_| PyValueError::new_err("each candidate must be a dict"))?;
        out.push(Candidate {
            left_seq: get_required_string(&dict, "left_seq")?,
            right_seq: get_required_string(&dict, "right_seq")?,
            product_size: get_required_i64(&dict, "product_size")?,
            left_len: get_required_i64(&dict, "left_len")?,
            right_len: get_required_i64(&dict, "right_len")?,
            tm_left: get_optional_f64(&dict, "tm_left")?,
            tm_right: get_optional_f64(&dict, "tm_right")?,
            gc_left: get_optional_f64(&dict, "gc_left")?,
            gc_right: get_optional_f64(&dict, "gc_right")?,
            pair_penalty: get_optional_f64(&dict, "pair_penalty")?,
            left_start: get_required_i64(&dict, "left_start")?,
        });
    }
    Ok(out)
}

fn has_3p_gc_run(seq: &str, run_len: usize) -> bool {
    let mut count = 0usize;
    for base in seq.chars().rev() {
        if base == 'G' || base == 'C' {
            count += 1;
            if count >= run_len {
                return true;
            }
        } else {
            break;
        }
    }
    false
}

fn has_poly_run(seq: &str, run_len: usize) -> bool {
    let mut chars = seq.chars();
    let Some(mut last) = chars.next() else {
        return false;
    };
    let mut count = 1usize;
    for ch in chars {
        if ch == last {
            count += 1;
            if count >= run_len {
                return true;
            }
        } else {
            last = ch;
            count = 1;
        }
    }
    false
}

fn junction_flags(
    boundaries: &[i64],
    product_start: i64,
    left_len: i64,
    right_len: i64,
    product_size: i64,
) -> (bool, String) {
    if boundaries.is_empty() {
        return (false, "no_exon_info".to_string());
    }

    let left_0 = product_start;
    let left_1 = product_start + left_len - 1;
    let right_3p = product_start + product_size - 1;

    let spans_left = boundaries.iter().any(|&b| left_0 <= b && b < left_1);
    let spans_amplicon = boundaries
        .iter()
        .any(|&b| product_start <= b && b < product_start + product_size - 1);

    if spans_left {
        return (true, "left_primer_spans_junction".to_string());
    }
    if spans_amplicon {
        return (true, "amplicon_spans_junction".to_string());
    }

    let right_start_guess = right_3p - (right_len - 1);
    let spans_right = boundaries
        .iter()
        .any(|&b| right_start_guess <= b && b < right_3p);
    if spans_right {
        return (true, "right_primer_spans_junction (approx)".to_string());
    }
    (false, "no_junction".to_string())
}

fn probe_preference(product_size: i64, probe_start: i64, probe_end: i64, dist_l: i64, dist_r: i64) -> f64 {
    let min_dist = dist_l.min(dist_r) as f64;
    let imbalance = (dist_l - dist_r).abs() as f64;
    let probe_mid = (probe_start + probe_end) as f64 / 2.0;
    let center = (product_size - 1) as f64 / 2.0;
    let center_offset = (probe_mid - center).abs();
    (min_dist * 0.2) - (imbalance * 0.5) - (center_offset * 0.2)
}

fn score_candidate(
    pair_penalty: Option<f64>,
    junction_spanning: bool,
    dist_l: i64,
    dist_r: i64,
    probe_mid: f64,
    product_size: i64,
) -> f64 {
    let mut score = 0.0;
    if let Some(penalty) = pair_penalty {
        score += (10.0 - penalty).max(0.0);
    }
    if junction_spanning {
        score += 5.0;
    }

    let min_dist = dist_l.min(dist_r) as f64;
    let imbalance = (dist_l - dist_r).abs() as f64;
    let center = (product_size - 1) as f64 / 2.0;
    let center_offset = (probe_mid - center).abs();

    score += min_dist.min(30.0) * 0.05;
    score -= imbalance.min(60.0) * 0.02;
    score -= center_offset.min(60.0) * 0.01;
    score
}

fn evaluate_candidate(
    cand: &Candidate,
    template: &str,
    boundaries: &[i64],
    min_probe_offset_bp: i64,
    ac: &AhoCorasick,
    probes: &[Probe],
    pattern_meta: &[PatternMeta],
) -> Option<RankedResult> {
    if cand.product_size <= 0 || cand.left_len <= 0 || cand.right_len <= 0 || cand.left_start < 0 {
        return None;
    }

    let product_start = cand.left_start;
    let product_end = product_start + cand.product_size;
    if product_end <= product_start {
        return None;
    }

    let ps = usize::try_from(product_start).ok()?;
    let pe = usize::try_from(product_end).ok()?;
    if pe > template.len() {
        return None;
    }

    let amplicon = &template[ps..pe];
    if amplicon.len() != cand.product_size as usize {
        return None;
    }

    if has_3p_gc_run(&cand.left_seq, 3) || has_3p_gc_run(&cand.right_seq, 3) {
        return None;
    }
    if has_poly_run(&cand.left_seq, 4) || has_poly_run(&cand.right_seq, 4) {
        return None;
    }

    let hits = collect_matches_from_normalized(ac, probes, pattern_meta, amplicon);
    if hits.is_empty() {
        return None;
    }

    let (junction_spanning, detail) = junction_flags(
        boundaries,
        product_start,
        cand.left_len,
        cand.right_len,
        cand.product_size,
    );

    let mut best_hit: Option<(f64, ProbeMatch, i64, i64)> = None;
    for hit in hits {
        // Amplicon-relative coordinates (0-based, inclusive)
        let left_3p = cand.left_len - 1;
        let right_3p = cand.product_size - 1;
        let right_primer_start = cand.product_size - cand.right_len;

        // Always forbid primer/probe overlap even when min_probe_offset_bp=0.
        if (hit.start as i64) <= left_3p {
            continue;
        }
        if (hit.end as i64) >= right_primer_start {
            continue;
        }

        // Distance in bp between primer 3' end and nearest probe base (0 when adjacent)
        let dist_l = hit.start as i64 - left_3p - 1;
        let dist_r = right_3p - hit.end as i64 - 1;
        if dist_l < min_probe_offset_bp || dist_r < min_probe_offset_bp {
            continue;
        }
        let pref = probe_preference(
            cand.product_size,
            hit.start as i64,
            hit.end as i64,
            dist_l,
            dist_r,
        );
        match &best_hit {
            Some((best_pref, _, _, _)) if pref <= *best_pref => {}
            _ => best_hit = Some((pref, hit, dist_l, dist_r)),
        }
    }

    let (_, hit, dist_l, dist_r) = best_hit?;

    let probe_mid = (hit.start + hit.end) as f64 / 2.0;
    let score = score_candidate(
        cand.pair_penalty,
        junction_spanning,
        dist_l,
        dist_r,
        probe_mid,
        cand.product_size,
    );

    let left_start = cand.left_start;
    let left_end = cand.left_start + cand.left_len - 1;
    let right_end = product_start + cand.product_size - 1;
    let right_start = right_end - cand.right_len + 1;

    Some(RankedResult {
        left_seq: cand.left_seq.clone(),
        right_seq: cand.right_seq.clone(),
        product_size: cand.product_size,
        tm_left: cand.tm_left,
        tm_right: cand.tm_right,
        gc_left: cand.gc_left,
        gc_right: cand.gc_right,
        pair_penalty: cand.pair_penalty,
        upl_probe_id: hit.probe_id,
        upl_probe_seq: hit.probe_seq,
        upl_probe_strand: hit.strand,
        probe_start_in_amplicon: hit.start as i64,
        probe_end_in_amplicon: hit.end as i64,
        dist_left_3p_to_probe: dist_l,
        dist_probe_to_right_3p: dist_r,
        left_start,
        left_end,
        right_start,
        right_end,
        left_len: cand.left_len,
        right_len: cand.right_len,
        amplicon_start: product_start,
        amplicon_end: product_end - 1,
        junction_spanning,
        junction_spanning_detail: detail,
        score,
    })
}

fn sort_ranked(a: &RankedResult, b: &RankedResult) -> std::cmp::Ordering {
    b.score
        .total_cmp(&a.score)
        .then_with(|| a.pair_penalty.is_none().cmp(&b.pair_penalty.is_none()))
        .then_with(|| {
            let ap = a.pair_penalty.unwrap_or(1e9);
            let bp = b.pair_penalty.unwrap_or(1e9);
            ap.total_cmp(&bp)
        })
}

fn rank_candidates_impl(
    template: &str,
    probes: &[Probe],
    candidates: &[Candidate],
    boundaries: &[i64],
    min_probe_offset_bp: i64,
) -> Vec<RankedResult> {
    if probes.is_empty() || candidates.is_empty() {
        return Vec::new();
    }

    let Ok((ac, pattern_meta)) = build_matcher(probes) else {
        return Vec::new();
    };
    let template = normalize_seq(template);

    let ranked: Vec<RankedResult> = candidates
        .par_iter()
        .filter_map(|cand| {
            evaluate_candidate(
                cand,
                &template,
                boundaries,
                min_probe_offset_bp,
                &ac,
                probes,
                &pattern_meta,
            )
        })
        .collect();

    let mut best: HashMap<(String, String, i64, i64, i64), RankedResult> = HashMap::new();
    for p in ranked {
        let key = (
            p.left_seq.clone(),
            p.right_seq.clone(),
            p.left_start,
            p.right_start,
            p.product_size,
        );
        let Some(prev) = best.get(&key) else {
            best.insert(key, p);
            continue;
        };

        if p.score > prev.score {
            best.insert(key, p);
            continue;
        }
        if p.score == prev.score {
            let p_min = p.dist_left_3p_to_probe.min(p.dist_probe_to_right_3p);
            let prev_min = prev
                .dist_left_3p_to_probe
                .min(prev.dist_probe_to_right_3p);
            if p_min > prev_min {
                best.insert(key, p);
            }
        }
    }

    let mut out: Vec<RankedResult> = best.into_values().collect();
    out.sort_by(sort_ranked);
    out
}

fn ranked_to_pydict(py: Python<'_>, r: &RankedResult) -> PyResult<Py<PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("left_seq", &r.left_seq)?;
    d.set_item("right_seq", &r.right_seq)?;
    d.set_item("product_size", r.product_size)?;
    d.set_item("tm_left", r.tm_left)?;
    d.set_item("tm_right", r.tm_right)?;
    d.set_item("gc_left", r.gc_left)?;
    d.set_item("gc_right", r.gc_right)?;
    d.set_item("pair_penalty", r.pair_penalty)?;
    d.set_item("upl_probe_id", &r.upl_probe_id)?;
    d.set_item("upl_probe_seq", &r.upl_probe_seq)?;
    d.set_item("upl_probe_strand", &r.upl_probe_strand)?;
    d.set_item("probe_start_in_amplicon", r.probe_start_in_amplicon)?;
    d.set_item("probe_end_in_amplicon", r.probe_end_in_amplicon)?;
    d.set_item("dist_left_3p_to_probe", r.dist_left_3p_to_probe)?;
    d.set_item("dist_probe_to_right_3p", r.dist_probe_to_right_3p)?;
    d.set_item("left_start", r.left_start)?;
    d.set_item("left_end", r.left_end)?;
    d.set_item("right_start", r.right_start)?;
    d.set_item("right_end", r.right_end)?;
    d.set_item("left_len", r.left_len)?;
    d.set_item("right_len", r.right_len)?;
    d.set_item("amplicon_start", r.amplicon_start)?;
    d.set_item("amplicon_end", r.amplicon_end)?;
    d.set_item("junction_spanning", r.junction_spanning)?;
    d.set_item("junction_spanning_detail", &r.junction_spanning_detail)?;
    d.set_item("score", r.score)?;
    Ok(d.unbind())
}

#[pymethods]
impl UplFinder {
    #[new]
    fn new(probes: &Bound<'_, PyAny>) -> PyResult<Self> {
        let parsed = parse_probes_from_py(probes)?;
        if parsed.is_empty() {
            return Ok(UplFinder {
                ac: None,
                probes: parsed,
                pattern_meta: Vec::new(),
            });
        }

        let (ac, pattern_meta) = build_matcher(&parsed).map_err(PyValueError::new_err)?;
        Ok(UplFinder {
            ac: Some(ac),
            probes: parsed,
            pattern_meta,
        })
    }

    fn find(&self, sequence: &str) -> Vec<ProbeMatch> {
        let Some(ac) = &self.ac else {
            return Vec::new();
        };
        let seq = normalize_seq(sequence);
        collect_matches_from_normalized(ac, &self.probes, &self.pattern_meta, &seq)
    }
}

#[pyfunction]
fn find_upl_matches_from_file(sequence: &str, probes_path: &str) -> PyResult<Vec<ProbeMatch>> {
    let probes = load_upl_probes(probes_path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(find_upl_matches(sequence, &probes))
}

#[pyfunction]
fn find_upl_matches_from_pickle(sequence: &str, probes_path: &str) -> PyResult<Vec<ProbeMatch>> {
    find_upl_matches_from_file(sequence, probes_path)
}

#[pyfunction(signature = (template, probes, candidates, exon_boundaries=None, min_probe_offset_bp=None))]
fn rank_candidates(
    py: Python<'_>,
    template: &str,
    probes: &Bound<'_, PyAny>,
    candidates: &Bound<'_, PyAny>,
    exon_boundaries: Option<Vec<i64>>,
    min_probe_offset_bp: Option<i64>,
) -> PyResult<Vec<Py<PyDict>>> {
    let parsed_probes = parse_probes_from_py(probes)?;
    let parsed_candidates = parse_candidates_from_py(candidates)?;
    let boundaries = exon_boundaries.unwrap_or_default();

    let ranked = rank_candidates_impl(
        template,
        &parsed_probes,
        &parsed_candidates,
        &boundaries,
        min_probe_offset_bp.unwrap_or(0),
    );

    ranked.iter().map(|r| ranked_to_pydict(py, r)).collect()
}

#[pymodule]
fn upl_finder_rust(_py: Python<'_>, m: &pyo3::Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_upl_matches_from_pickle, m)?)?;
    m.add_function(wrap_pyfunction!(find_upl_matches_from_file, m)?)?;
    m.add_function(wrap_pyfunction!(rank_candidates, m)?)?;
    m.add_class::<UplFinder>()?;
    m.add_class::<ProbeMatch>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_works() {
        assert_eq!(revcomp("ACTGN"), "NCAGT");
    }

    #[test]
    fn find_matches_supports_strand_and_palindrome() {
        let probes = vec![
            Probe {
                probe_id: "p1".to_string(),
                probe_seq: "ACTG".to_string(),
            },
            Probe {
                probe_id: "pal".to_string(),
                probe_seq: "AGCT".to_string(),
            },
            Probe {
                probe_id: "ov".to_string(),
                probe_seq: "AAA".to_string(),
            },
        ];
        let hits = find_upl_matches("ACTGCAGTAAAAAGCT", &probes);

        assert!(hits.iter().any(|h| h.probe_id == "p1" && h.start == 0 && h.strand == "+"));
        assert!(hits.iter().any(|h| h.probe_id == "p1" && h.start == 4 && h.strand == "-"));
        assert!(hits
            .iter()
            .filter(|h| h.probe_id == "pal")
            .all(|h| h.strand == "+"));
        assert!(hits
            .iter()
            .filter(|h| h.probe_id == "ov")
            .count()
            >= 3);
    }

    #[test]
    fn rank_candidates_returns_probe_and_score_fields() {
        let probes = vec![Probe {
            probe_id: "p1".to_string(),
            probe_seq: "CCCG".to_string(),
        }];
        let candidates = vec![Candidate {
            left_seq: "ATATATAT".to_string(),
            right_seq: "TATATATA".to_string(),
            product_size: 24,
            left_len: 4,
            right_len: 4,
            tm_left: Some(60.0),
            tm_right: Some(60.0),
            gc_left: Some(50.0),
            gc_right: Some(50.0),
            pair_penalty: Some(1.5),
            left_start: 0,
        }];

        let ranked = rank_candidates_impl(
            "TTTAAACCCGGGTTTAAACCCGGG",
            &probes,
            &candidates,
            &[],
            2,
        );

        assert_eq!(ranked.len(), 1);
        let first = &ranked[0];
        assert_eq!(first.upl_probe_id, "p1");
        assert!(first.upl_probe_strand == "+" || first.upl_probe_strand == "-");
        assert!(first.probe_start_in_amplicon == 6 || first.probe_start_in_amplicon == 8);
        assert!(first.score > 0.0);
    }

    #[test]
    fn rank_candidates_rejects_probe_overlapping_primers_even_when_offset_zero() {
        let probes = vec![Probe {
            probe_id: "p1".to_string(),
            probe_seq: "CCCC".to_string(),
        }];
        let candidates = vec![Candidate {
            left_seq: "ATATAT".to_string(),
            right_seq: "TATATA".to_string(),
            product_size: 20,
            left_len: 6,
            right_len: 6,
            tm_left: Some(60.0),
            tm_right: Some(60.0),
            gc_left: Some(50.0),
            gc_right: Some(50.0),
            pair_penalty: Some(1.0),
            left_start: 0,
        }];

        // Overlaps left primer region [0..5]
        let ranked_left = rank_candidates_impl("AACCCCAAAAAAAAAAAAAA", &probes, &candidates, &[], 0);
        assert_eq!(ranked_left.len(), 0);

        // Overlaps right primer region [14..19]
        let ranked_right = rank_candidates_impl("AAAAAAAAAAAAAAACCCCA", &probes, &candidates, &[], 0);
        assert_eq!(ranked_right.len(), 0);

        // Internal probe is allowed when offset=0 (adjacent is 0bp, overlap is rejected above)
        let ranked_ok = rank_candidates_impl("AAAAAAAACCCCAAAAAAAA", &probes, &candidates, &[], 0);
        assert_eq!(ranked_ok.len(), 1);
        assert_eq!(ranked_ok[0].probe_start_in_amplicon, 8);
        assert_eq!(ranked_ok[0].probe_end_in_amplicon, 11);
        assert_eq!(ranked_ok[0].dist_left_3p_to_probe, 2);
        assert_eq!(ranked_ok[0].dist_probe_to_right_3p, 7);
    }

    #[test]
    fn parse_probes_from_tsv_text_parses_and_normalizes() {
        let text = r#"
# comment
p1	 ac tg
p2	NNNN
"#;
        let probes = parse_probes_from_tsv_text(text).unwrap();
        assert_eq!(probes.len(), 2);
        assert_eq!(probes[0].probe_id, "p1");
        assert_eq!(probes[0].probe_seq, "ACTG");
        assert_eq!(probes[1].probe_id, "p2");
        assert_eq!(probes[1].probe_seq, "NNNN");
    }

    #[test]
    fn parse_probes_from_json_text_parses_and_normalizes() {
        let text = r#"{"p1":"ac tg","p2":"NNNN"}"#;
        let mut probes = parse_probes_from_json_text(text).unwrap();
        probes.sort_by(|a, b| a.probe_id.cmp(&b.probe_id));
        assert_eq!(probes.len(), 2);
        assert_eq!(probes[0].probe_id, "p1");
        assert_eq!(probes[0].probe_seq, "ACTG");
        assert_eq!(probes[1].probe_id, "p2");
        assert_eq!(probes[1].probe_seq, "NNNN");
    }
}
