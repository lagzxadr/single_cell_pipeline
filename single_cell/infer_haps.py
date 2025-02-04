'''
Created on Aug 29, 2018

@author: dgrewal
'''

import os
import sys

import pypeliner.managed as mgd
from single_cell.utils import inpututils

import pypeliner


def infer_haps(
        bam_file,
        haplotypes_filename,
        allele_counts_filename,
        config,
        from_tumour=False,
):
    baseimage = {'docker_image': config['docker']['single_cell_pipeline']}

    remixt_config = config.get('extract_seqdata', {})
    remixt_ref_data_dir = config['ref_data_dir']

    chromosomes = config['chromosomes']

    ctx = dict(mem_retry_increment=2, disk_retry_increment=50, ncpus=1, **baseimage)
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    if isinstance(bam_file, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('cell_id'),
            value=list(bam_file.keys()),
        )

        # dont parallelize over chromosomes for per cell bams
        workflow.subworkflow(
            name="extract_seqdata",
            axes=('cell_id',),
            func='single_cell.workflows.extract_seqdata.create_extract_seqdata_workflow',
            args=(
                mgd.InputFile(
                    'bam_markdups',
                    'cell_id',
                    fnames=bam_file,
                    extensions=['.bai']
                ),
                mgd.TempOutputFile('seqdata_cell.h5', 'cell_id'),
                config.get('extract_seqdata', {}),
                config['ref_data_dir'],
                config,
            )
        )
        workflow.transform(
            name='merge_all_seqdata',
            func="single_cell.workflows.titan.tasks.merge_overlapping_seqdata",
            args=(
                mgd.TempOutputFile('seqdata_file.h5'),
                mgd.TempInputFile("seqdata_cell.h5", "cell_id"),
                config["chromosomes"]
            ),
        )

    else:
        # if its a single bam, then its probably whole genome
        # so parallelize over chromosomes
        workflow.subworkflow(
            name='extract_seqdata',
            func='remixt.workflow.create_extract_seqdata_workflow',
            ctx={'disk': 150},
            args=(
                mgd.InputFile(bam_file, extensions=['.bai']),
                mgd.TempOutputFile('seqdata_file.h5'),
                remixt_config,
                remixt_ref_data_dir,
            ),
        )

    workflow.setobj(
        obj=mgd.OutputChunks('chromosome'),
        value=chromosomes,
    )

    if from_tumour:
        func = 'remixt.analysis.haplotype.infer_snp_genotype_from_tumour'
    else:
        func = 'remixt.analysis.haplotype.infer_snp_genotype_from_normal'

    workflow.transform(
        name='infer_snp_genotype',
        axes=('chromosome',),
        ctx=dict(mem=16, **ctx),
        func=func,
        args=(
            mgd.TempOutputFile('snp_genotype.tsv', 'chromosome'),
            mgd.TempInputFile('seqdata_file.h5'),
            mgd.InputInstance('chromosome'),
            config,
        ),
    )

    workflow.transform(
        name='infer_haps',
        axes=('chromosome',),
        ctx=dict(mem=16, **ctx),
        func='remixt.analysis.haplotype.infer_haps',
        args=(
            mgd.TempOutputFile('haplotypes.tsv', 'chromosome'),
            mgd.TempInputFile('snp_genotype.tsv', 'chromosome'),
            mgd.InputInstance('chromosome'),
            mgd.TempSpace('haplotyping', 'chromosome'),
            remixt_config,
            remixt_ref_data_dir,
        ),
    )

    workflow.transform(
        name='merge_haps',
        ctx=dict(mem=16, **ctx),
        func='remixt.utils.merge_tables',
        args=(
            mgd.OutputFile(haplotypes_filename),
            mgd.TempInputFile('haplotypes.tsv', 'chromosome'),
        )
    )

    workflow.transform(
        name='create_segments',
        ctx=dict(mem=16, **ctx),
        func='remixt.analysis.segment.create_segments',
        args=(
            mgd.TempOutputFile('segments.tsv'),
            remixt_config,
            config['ref_data_dir'],
        ),
    )

    workflow.transform(
        name='haplotype_allele_readcount',
        ctx=dict(mem=16, **ctx),
        func='remixt.analysis.readcount.haplotype_allele_readcount',
        args=(
            mgd.OutputFile(allele_counts_filename),
            mgd.TempInputFile('segments.tsv'),
            mgd.TempInputFile('seqdata_file.h5'),
            mgd.InputFile(haplotypes_filename),
            remixt_config
        ),
    )

    return workflow


def extract_allele_readcounts(
        haplotypes_filename,
        cell_bams,
        allele_counts_filename,
        config,
):
    baseimage = {'docker_image': config['docker']['single_cell_pipeline']}

    remixt_config = config.get('extract_seqdata', {})
    remixt_ref_data_dir = config['ref_data_dir']

    workflow = pypeliner.workflow.Workflow(ctx=baseimage)

    workflow.set_filenames('cell.bam', 'cell_id', fnames=cell_bams)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(cell_bams.keys()),
    )

    workflow.subworkflow(
        name='create_chromosome_seqdata',
        axes=('cell_id',),
        func='single_cell.workflows.extract_seqdata.create_extract_seqdata_workflow',
        args=(
            mgd.InputFile('cell.bam', 'cell_id', extensions=['.bai']),
            mgd.TempOutputFile('seqdata.h5', 'cell_id', axes_origin=[]),
            config.get('extract_seqdata', {}),
            config['ref_data_dir'],
            config,
        ),
    )

    # TODO Segments with bin width from single cell
    workflow.transform(
        name='create_segments',
        func='remixt.analysis.segment.create_segments',
        args=(
            mgd.TempOutputFile('segments.tsv'),
            remixt_config,
            remixt_ref_data_dir,
        ),
    )

    workflow.transform(
        name='haplotype_allele_readcount',
        axes=('cell_id',),
        ctx={'mem': 16},
        func='remixt.analysis.readcount.haplotype_allele_readcount',
        args=(
            mgd.TempOutputFile('allele_counts.tsv', 'cell_id', axes_origin=[]),
            mgd.TempInputFile('segments.tsv'),
            mgd.TempInputFile('seqdata.h5', 'cell_id'),
            mgd.InputFile(haplotypes_filename),
            remixt_config,
        ),
    )

    workflow.transform(
        name='merge_allele_readcount',
        ctx={'mem': 16},
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            mgd.TempInputFile('allele_counts.tsv', 'cell_id'),
            mgd.OutputFile(allele_counts_filename, extensions=['.yaml']),
        ),
        kwargs={
            'key_column': 'cell_id',
            'write_header': True
        },
    )

    return workflow


def infer_haps_workflow(args):
    config = inpututils.load_config(args)
    config = config['infer_haps']
    baseimage = config['docker']['single_cell_pipeline']

    ctx = dict(mem_retry_increment=2, disk_retry_increment=50, ncpus=1, docker_image=baseimage)
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    haplotypes_filename = os.path.join(args["out_dir"], "haplotypes.tsv")
    allele_counts_filename = os.path.join(args["out_dir"], "allele_counts.tsv")

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    normal_data, tumour_cells = inpututils.load_haps_input(args['input_yaml'])

    if isinstance(normal_data, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('normal_cell_id'),
            value=list(normal_data.keys()),
        )
        bam_file = mgd.InputFile('normal.bam', 'normal_cell_id', fnames=normal_data, extensions=['.bai'])
    else:
        bam_file = mgd.InputFile(normal_data, extensions=['.bai'])

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=list(tumour_cells.keys()),
    )

    workflow.subworkflow(
        name='infer_haps',
        func=infer_haps,
        args=(
            bam_file,
            mgd.OutputFile(haplotypes_filename),
            mgd.TempOutputFile("allele_counts.csv"),
            config,
        ),
    )

    workflow.subworkflow(
        name='extract_allele_readcounts',
        func=extract_allele_readcounts,
        args=(
            mgd.InputFile(haplotypes_filename),
            mgd.InputFile('tumour_cells.bam', 'tumour_cell_id', extensions=['.bai'],
                          axes_origin=[], fnames=tumour_cells),
            mgd.OutputFile(allele_counts_filename),
            config,
        ),
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            [haplotypes_filename, allele_counts_filename],
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'infer_haps'}
        }
    )

    return workflow


def infer_haps_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = infer_haps_workflow(args)

    pyp.run(workflow)
