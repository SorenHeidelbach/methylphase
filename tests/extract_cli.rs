use std::{
    ffi::OsString,
    fs,
    fs::File,
    io::Cursor,
    path::{Path, PathBuf},
};

use anyhow::Result;
use methylphase::{
    cli::{Cli, Command, ContigSelection},
    run,
};
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as SamWrite};
use tempfile::tempdir;

#[test]
fn extract_command_produces_expected_reports() -> Result<()> {
    let tmp = tempdir()?;
    let bam_path = tmp.path().join("modcalls.bam");
    write_test_bam(&bam_path)?;

    let output_dir = tmp.path().join("outputs_explicit");

    let cli = Cli {
        command: Command::Extract {
            bam: bam_path,
            motifs: vec!["GATC_6mA_1".to_string()],
            motif_file: None,
            motif_summary_tsv: None,
            fastq_dir: None,
            sequence_fallback: None,
            sequence_index: None,
            output_dir: output_dir.clone(),
            contig_args: ContigSelection::default(),
        },
    };

    run(cli)?;

    let per_contents = fs::read_to_string(output_dir.join("per_read.tsv"))?;
    let aggregate_contents = fs::read_to_string(output_dir.join("aggregate.tsv"))?;

    let expected_per = "\
read_id\tcontig\tstrand\tmotif\tmotif_start\tmotif_position\tread_position\tref_position\tprobability\tmethylated\tmod_label
read0\tctg\t+\tGATC_6mA_1\t0\t1\t1\t2\t0.7843\t1\t6mA
";

    let expected_aggregate = "\
read_id\tmotif\tcall_count\tmean_probability\tquantile_0\tquantile_10\tquantile_20\tquantile_30\tquantile_40\tquantile_50\tquantile_60\tquantile_70\tquantile_80\tquantile_90\tquantile_100
read0\tGATC_6mA_1\t1\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843
";

    assert_eq!(per_contents, expected_per);
    assert_eq!(aggregate_contents, expected_aggregate);

    Ok(())
}

#[test]
fn extract_command_can_write_to_output_dir_defaults() -> Result<()> {
    let tmp = tempdir()?;
    let bam_path = tmp.path().join("modcalls.bam");
    write_test_bam(&bam_path)?;
    let out_dir = tmp.path().join("outputs");

    let cli = Cli {
        command: Command::Extract {
            bam: bam_path,
            motifs: vec!["GATC_6mA_1".to_string()],
            motif_file: None,
            motif_summary_tsv: None,
            fastq_dir: None,
            sequence_fallback: None,
            sequence_index: None,
            output_dir: out_dir.clone(),
            contig_args: ContigSelection::default(),
        },
    };

    run(cli)?;

    let per_contents = fs::read_to_string(out_dir.join("per_read.tsv"))?;
    let aggregate_contents = fs::read_to_string(out_dir.join("aggregate.tsv"))?;

    let expected_per = "\
read_id\tcontig\tstrand\tmotif\tmotif_start\tmotif_position\tread_position\tref_position\tprobability\tmethylated\tmod_label
read0\tctg\t+\tGATC_6mA_1\t0\t1\t1\t2\t0.7843\t1\t6mA
";

    let expected_aggregate = "\
read_id\tmotif\tcall_count\tmean_probability\tquantile_0\tquantile_10\tquantile_20\tquantile_30\tquantile_40\tquantile_50\tquantile_60\tquantile_70\tquantile_80\tquantile_90\tquantile_100
read0\tGATC_6mA_1\t1\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843
";

    assert_eq!(per_contents, expected_per);
    assert_eq!(aggregate_contents, expected_aggregate);
    Ok(())
}

#[test]
fn split_reads_command_clusters_and_writes_fastqs() -> Result<()> {
    let tmp = tempdir()?;
    let bam_path = tmp.path().join("modcalls.bam");
    write_test_bam(&bam_path)?;
    let out_dir = tmp.path().join("split_outputs");

    let cli = Cli {
        command: Command::SplitReads {
            bam: bam_path,
            motifs: vec!["GATC_6mA_1".to_string()],
            motif_file: None,
            sequence_fallback: None,
            sequence_index: None,
            output_dir: out_dir.clone(),
            min_cluster_size: 1,
            min_samples: Some(1),
            emit_fastq: true,
            threads: 1,
            contig_args: ContigSelection::default(),
        },
    };

    run(cli)?;

    let clustering = fs::read_to_string(out_dir.join("read_clustering.tsv"))?;
    assert!(clustering.contains("read0"));
    let raw_clustering = fs::read_to_string(out_dir.join("read_clustering_raw.tsv"))?;
    assert!(raw_clustering.contains("read0"));
    let summary = fs::read_to_string(out_dir.join("cluster_summary.tsv"))?;
    assert!(summary.contains("cluster_id"));
    let clusters_dir = out_dir.join("clusters");
    let noise_fastq = clusters_dir.join("cluster_noise.fastq");
    assert!(noise_fastq.exists());
    Ok(())
}

#[test]
fn split_reads_supports_multiple_threads() -> Result<()> {
    let tmp = tempdir()?;
    let bam_path = tmp.path().join("modcalls.bam");
    write_test_bam(&bam_path)?;
    let out_dir = tmp.path().join("split_outputs_threads");

    let cli = Cli {
        command: Command::SplitReads {
            bam: bam_path,
            motifs: vec!["GATC_6mA_1".to_string()],
            motif_file: None,
            sequence_fallback: None,
            sequence_index: None,
            output_dir: out_dir.clone(),
            min_cluster_size: 1,
            min_samples: Some(1),
            emit_fastq: true,
            threads: 2,
            contig_args: ContigSelection::default(),
        },
    };

    run(cli)?;
    assert!(out_dir.join("read_clustering.tsv").exists());
    assert!(out_dir.join("read_clustering_raw.tsv").exists());
    Ok(())
}

#[test]
fn split_reads_imputes_missing_data() -> Result<()> {
    let tmp = tempdir()?;
    let bam_path = tmp.path().join("modcalls.bam");
    write_test_bam(&bam_path)?;
    let out_dir = tmp.path().join("split_outputs_missing");

    let cli = Cli {
        command: Command::SplitReads {
            bam: bam_path,
            motifs: vec!["GATC_6mA_1".to_string(), "CCWGG_5mC_1".to_string()],
            motif_file: None,
            sequence_fallback: None,
            sequence_index: None,
            output_dir: out_dir.clone(),
            min_cluster_size: 1,
            min_samples: Some(1),
            emit_fastq: true,
            threads: 1,
            contig_args: ContigSelection::default(),
        },
    };

    run(cli)?;
    let clustering = fs::read_to_string(out_dir.join("read_clustering.tsv"))?;
    assert!(clustering.contains("GATC_6mA_1"));
    assert!(!clustering.contains("CCWGG_5mC_1"));
    let raw_clustering = fs::read_to_string(out_dir.join("read_clustering_raw.tsv"))?;
    assert!(!raw_clustering.contains("NA"));
    Ok(())
}

#[test]
fn extract_handles_ambiguous_motif() -> Result<()> {
    let tmp = tempdir()?;
    let bam_path = tmp.path().join("degenerate.bam");
    write_degenerate_bam(&bam_path)?;

    let output_dir = tmp.path().join("out");
    let cli = Cli {
        command: Command::Extract {
            bam: bam_path,
            motifs: vec!["AAANNGTG_6mA_1".to_string()],
            motif_file: None,
            motif_summary_tsv: None,
            fastq_dir: None,
            sequence_fallback: None,
            sequence_index: None,
            output_dir: output_dir.clone(),
            contig_args: ContigSelection::default(),
        },
    };

    run(cli)?;

    let per_path = output_dir.join("per_read.tsv");
    let agg_path = output_dir.join("aggregate.tsv");
    let per_contents = fs::read_to_string(per_path)?;
    let agg_contents = fs::read_to_string(agg_path)?;

    let expected_per = "\
read_id\tcontig\tstrand\tmotif\tmotif_start\tmotif_position\tread_position\tref_position\tprobability\tmethylated\tmod_label
read0\tctg\t+\tAAANNGTG_6mA_1\t0\t1\t1\t2\t0.7843\t1\t6mA
";
    let expected_agg = "\
read_id\tmotif\tcall_count\tmean_probability\tquantile_0\tquantile_10\tquantile_20\tquantile_30\tquantile_40\tquantile_50\tquantile_60\tquantile_70\tquantile_80\tquantile_90\tquantile_100
read0\tAAANNGTG_6mA_1\t1\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843
";

    assert_eq!(per_contents, expected_per);
    assert_eq!(agg_contents, expected_agg);
    Ok(())
}

#[test]
fn extract_handles_reverse_read() -> Result<()> {
    let tmp = tempdir()?;
    let bam_path = tmp.path().join("reverse.bam");
    write_reverse_bam(&bam_path)?;
    let out_dir = tmp.path().join("reverse_out");

    let cli = Cli {
        command: Command::Extract {
            bam: bam_path,
            motifs: vec!["GATC_6mA_1".to_string()],
            motif_file: None,
            motif_summary_tsv: None,
            fastq_dir: None,
            sequence_fallback: None,
            sequence_index: None,
            output_dir: out_dir.clone(),
            contig_args: ContigSelection::default(),
        },
    };

    run(cli)?;
    let per_contents = fs::read_to_string(out_dir.join("per_read.tsv"))?;
    assert!(per_contents.contains("read_rev"));
    Ok(())
}

#[test]
fn extract_filters_motifs_per_bin() -> Result<()> {
    let tmp = tempdir()?;
    let bam_path = tmp.path().join("modcalls.bam");
    write_test_bam(&bam_path)?;

    let bins_path = tmp.path().join("bins.tsv");
    fs::write(&bins_path, "ctg\tbinA\nctg\tbinB\n")?;

    let motif_path = tmp.path().join("motifs.tsv");
    fs::write(
        &motif_path,
        "id\tmotif\tmod_type\tmod_position\nbinA\tGATC\t6mA\t1\nbinB\tCCWGG\tm\t0\n",
    )?;

    let output_dir = tmp.path().join("bin_outputs");
    let cli = Cli {
        command: Command::Extract {
            bam: bam_path,
            motifs: Vec::new(),
            motif_file: Some(motif_path),
            motif_summary_tsv: None,
            fastq_dir: None,
            sequence_fallback: None,
            sequence_index: None,
            output_dir: output_dir.clone(),
            contig_args: ContigSelection {
                contigs: Vec::new(),
                contig_bins: Some(bins_path),
            },
        },
    };

    run(cli)?;

    let bin_a = fs::read_to_string(output_dir.join("binA").join("per_read.tsv"))?;
    assert!(bin_a.contains("GATC_6mA_1"));

    let bin_b = fs::read_to_string(output_dir.join("binB").join("per_read.tsv"))?;
    assert!(!bin_b.contains("GATC_6mA_1"));
    Ok(())
}

#[test]
fn vcf_command_writes_records() -> Result<()> {
    let tmp = tempdir()?;
    let bam_path = tmp.path().join("modcalls.bam");
    write_test_bam(&bam_path)?;
    let vcf_path = tmp.path().join("calls.vcf");

    let cli = Cli {
        command: Command::Vcf {
            bam: bam_path,
            motifs: vec!["GATC_6mA_1".to_string()],
            motif_file: None,
            sequence_fallback: None,
            sequence_index: None,
            output: vcf_path.clone(),
            sample_name: Some("sample1".to_string()),
            methylation_threshold: 0.5,
            plain_output: true,
            contig_args: ContigSelection::default(),
        },
    };

    run(cli)?;

    let contents = fs::read_to_string(&vcf_path)?;
    assert!(contents.starts_with("##fileformat=VCFv4.3"));
    assert!(contents.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"));
    assert!(contents.contains("ctg\t2\t.\tA\t<6mA>\t.\tPASS"));
    assert!(contents.contains("DP:METH:UNMETH:FRAC:MEANP"));
    assert!(contents.contains("sample1"));
    Ok(())
}

#[test]
fn impute_bam_rewrites_sequence_and_quality() -> Result<()> {
    let tmp = tempdir()?;
    let bam_path = tmp.path().join("modcalls.bam");
    write_test_bam(&bam_path)?;
    let output_path = tmp.path().join("imputed.bam");

    let cli = Cli {
        command: Command::ImputeBam {
            bam: bam_path.clone(),
            output: output_path.clone(),
            methylation_threshold: 0.5,
            summary: None,
            motifs: Vec::new(),
            motif_file: None,
            impute_all: true,
            contig_args: ContigSelection::default(),
        },
    };

    run(cli)?;

    let mut reader = bam::io::reader::Builder::default().build_from_path(&output_path)?;
    reader.read_header()?;
    let mut record = bam::Record::default();
    assert!(reader.read_record(&mut record)? > 0);

    let seq: String = record
        .sequence()
        .iter()
        .map(|base| char::from(base))
        .collect();
    assert_eq!(seq, "GCTC");

    let quality_scores = record.quality_scores();
    assert!(!quality_scores.is_empty());
    let q_char = (quality_scores.as_ref()[1] + 33) as char;
    let expected = phred_char_from_probability(200.0 / 255.0);
    assert_eq!(q_char, expected);
    Ok(())
}

fn write_test_bam(path: &Path) -> Result<()> {
    const SAM_TEMPLATE: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ctg\tLN:50
read0\t0\tctg\t1\t60\t4M\t*\t0\t0\tGATC\tFFFF\tMM:Z:A+a,0;\tML:B:C,200
";

    let mut reader = sam::io::Reader::new(Cursor::new(SAM_TEMPLATE.as_bytes().to_vec()));
    let header = reader.read_header()?;

    let file = File::create(path)?;
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header)?;

    let mut record = sam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        writer.write_alignment_record(&header, &record)?;
    }

    writer.try_finish()?;

    let index = bam::fs::index(path)?;

    let mut bai_path = OsString::from(path.as_os_str());
    bai_path.push(".bai");
    bam::bai::fs::write(PathBuf::from(bai_path), &index)?;

    Ok(())
}

fn write_degenerate_bam(path: &Path) -> Result<()> {
    const SAM_TEMPLATE: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ctg\tLN:50
read0\t0\tctg\t1\t60\t8M\t*\t0\t0\tAAACCGTG\tFFFFFFFF\tMM:Z:A+a,1;\tML:B:C,200
";

    let mut reader = sam::io::Reader::new(Cursor::new(SAM_TEMPLATE.as_bytes().to_vec()));
    let header = reader.read_header()?;

    let file = File::create(path)?;
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header)?;

    let mut record = sam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        writer.write_alignment_record(&header, &record)?;
    }

    writer.try_finish()?;
    let index = bam::fs::index(path)?;
    let mut bai_path = OsString::from(path.as_os_str());
    bai_path.push(".bai");
    bam::bai::fs::write(PathBuf::from(bai_path), &index)?;
    Ok(())
}

fn write_reverse_bam(path: &Path) -> Result<()> {
    const SAM_TEMPLATE: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ctg\tLN:50
read_rev\t16\tctg\t1\t60\t4M\t*\t0\t0\tGATC\tFFFF\tMM:Z:A+a,0;\tML:B:C,200
";

    let mut reader = sam::io::Reader::new(Cursor::new(SAM_TEMPLATE.as_bytes().to_vec()));
    let header = reader.read_header()?;

    let file = File::create(path)?;
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header)?;

    let mut record = sam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        writer.write_alignment_record(&header, &record)?;
    }

    writer.try_finish()?;
    let index = bam::fs::index(path)?;
    let mut bai_path = OsString::from(path.as_os_str());
    bai_path.push(".bai");
    bam::bai::fs::write(PathBuf::from(bai_path), &index)?;
    Ok(())
}

fn phred_char_from_probability(prob: f32) -> char {
    let q = if prob >= 1.0 {
        93.0
    } else {
        let complement = (1.0 - prob).max(1e-9);
        (-10.0 * complement.log10()).round().clamp(0.0, 93.0)
    };
    let phred = q as u8;
    (phred + 33) as char
}
