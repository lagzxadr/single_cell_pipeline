#Douglas Abrams
#10/14/19

import pypeliner
import pypeliner.managed as mgd

'''
Calls task functions to genotype variant calls

Runs svtyper genotyping on lumpy/destruct variant calls
for multiple libraries/samples in parallel

Arguments:
    reference (str) path to reference fastas
    input_bam (pypeliner.mgd) input bam file path
    lumpy_csv (pypeliner.mgd) lumpy calls file path
    destruct_csv (pypeliner.mgd) destruct calls file path
    AO ... AB (pypeliner.mgd) output paths for svtyper annotations
    config (dict) config information for pypeline

'''
def genotyping(
        reference,
        input_bam,
        lumpy_csv,
        destruct_csv,
        AO, AP, AS,
        ASC,DP,GQ,
        QA,QR,RO,
        RP,RS,SQ,
        GL,AB,
        config      
):
    img = config["docker"]["scp"]
    genotyping_ctx = {'docker_image': img}

    workflow = pypeliner.workflow.Workflow(ctx=genotyping_ctx)


    #lumpy csv -> vcf
    workflow.transform(
        name='lumpy_csv_to_vcf',
        func='single_cell.workflows.genotyping.tasks.varcalls_to_svtyper_input',
        args=(
            mgd.InputFile(lumpy_csv),
            mgd.TempOutputFile("lumpy_vcf"),
            mgd.TempSpace("lumpy_input"),
            "lumpy",
        )
    ) 

    #destruct csv -> vcf
    workflow.transform(
        name='destruct_csv_to_vcf',

        func='single_cell.workflows.genotyping.tasks.varcalls_to_svtyper_input',
        args=(
            mgd.InputFile(destruct_csv),
            mgd.TempOutputFile("destruct_vcf"),
            mgd.TempSpace("destruct_input"),
            "destruct",
        )
    )
    
    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(input_bam.keys()),
    )

    
    workflow.transform(
        name='genotype-lumpy',
        func='single_cell.workflows.genotyping.tasks.genotype',
        axes=('cell_id',),
        args=(
            mgd.InputFile("cell_specific.bam", "cell_id", fnames = input_bam, extensions=['.bai']),
            reference,
            mgd.TempInputFile("lumpy_vcf"),
            mgd.TempOutputFile("genotypedlumpy_vcf", "cell_id"),
            mgd.TempOutputFile("genotypedlumpy_csv", "cell_id"),
            mgd.TempSpace("temp", 'cell_id'),
            config["docker"]["svtyper"]
        )
    )

    workflow.transform(
        name='merge_cell_specific_csv-lumpy',
        func="single_cell.utils.vcfutils.merge_csvs",
        args=(
            mgd.TempInputFile("genotypedlumpy_csv", 'cell_id', axes_origin=[]),
            mgd.TempOutputFile("genotypedmerged_lumpy"),
            mgd.TempSpace("merged_lumpy")
        )
    )
    
    workflow.transform(
        name='genotype-destruct',
        func='single_cell.workflows.genotyping.tasks.genotype',
        axes=('cell_id',),
        args=(
            mgd.InputFile("cell_specific.bam", "cell_id", fnames = input_bam, extensions=['.bai']),
            reference,
            mgd.TempInputFile("destruct_vcf"),
            mgd.TempOutputFile("genotypeddestruct_vcf", "cell_id"),
            mgd.TempOutputFile("genotypeddestruct_csv", "cell_id"),
            mgd.TempSpace("tempdestruct", 'cell_id'),
            config["docker"]["svtyper"]
        )
    )

    workflow.transform(
        name='merge_cell_specific_csv-destruct',
        func="single_cell.utils.vcfutils.merge_csvs",
        args=(
            mgd.TempInputFile("genotypeddestruct_csv", 'cell_id', axes_origin=[]),
            mgd.TempOutputFile("genotypedmerged_destruct"),
            mgd.TempSpace("merged_destruct")
        )
    )
    

    workflow.transform(
        name='merge_lumpy_destruct_vcfs',
        func='single_cell.utils.vcfutils.merge_csvs',
        args=(
            [mgd.TempInputFile("genotypedmerged_lumpy"), 
                mgd.TempInputFile("genotypedmerged_destruct")],
            mgd.TempOutputFile("lumpy_destruct_merged.csv"),
            mgd.TempSpace("temp")
        )         
    )
    
    
    workflow.transform(
        name='write_svtyper_annotations',
        func='single_cell.workflows.genotyping.tasks.writeSvtyperAnnotations',
        args=(
            mgd.TempInputFile("lumpy_destruct_merged.csv"),
            [mgd.OutputFile(AO),
            mgd.OutputFile(AP),
            mgd.OutputFile(AS),
            mgd.OutputFile(ASC),
            mgd.OutputFile(DP),
            mgd.OutputFile(GQ),
            mgd.OutputFile(QA),
            mgd.OutputFile(QR),
            mgd.OutputFile(RO),
            mgd.OutputFile(RP),
            mgd.OutputFile(RS),
            mgd.OutputFile(SQ),
            mgd.OutputFile(GL),
            mgd.OutputFile(AB)],
            mgd.TempSpace("combined") 
        )
    )
    
    return workflow
    
