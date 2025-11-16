use std::{
    ffi::OsString,
    fs,
    fs::File,
    io::Cursor,
    path::{Path, PathBuf},
};

use anyhow::Result;
use methylation_phasing::{
    cli::{Cli, Command},
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

    let per_read = tmp.path().join("per.tsv");
    let aggregate = tmp.path().join("agg.tsv");

    let cli = Cli {
        command: Command::Extract {
            bam: bam_path,
            motifs: vec!["GATC_6mA_2".to_string()],
            motif_file: None,
            motif_summary_tsv: None,
            fastq_dir: None,
            per_read_tsv: per_read.clone(),
            aggregate_tsv: aggregate.clone(),
            contigs: Vec::new(),
        },
    };

    run(cli)?;

    let per_contents = fs::read_to_string(&per_read)?;
    let aggregate_contents = fs::read_to_string(&aggregate)?;

    let expected_per = "\
read_id\tcontig\tstrand\tmotif\tmotif_start\tmotif_position\tread_position\tref_position\tprobability\tmod_label
read0\tctg\t+\tGATC_6mA_2\t0\t2\t1\t2\t0.7843\t6mA
";

    let expected_aggregate = "\
read_id\tmotif\tcall_count\tmean_probability\tquantile_0\tquantile_10\tquantile_20\tquantile_30\tquantile_40\tquantile_50\tquantile_60\tquantile_70\tquantile_80\tquantile_90\tquantile_100
read0\tGATC_6mA_2\t1\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843\t0.7843
";

    assert_eq!(per_contents, expected_per);
    assert_eq!(aggregate_contents, expected_aggregate);

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
