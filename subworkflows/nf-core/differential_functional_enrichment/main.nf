
//
// Perform enrichment analysis
//
include { GPROFILER2_GOST          } from "../../../modules/nf-core/gprofiler2/gost/main.nf"
include { CUSTOM_TABULARTOGSEAGCT  } from '../../../modules/nf-core/custom/tabulartogseagct/main.nf'
include { CUSTOM_TABULARTOGSEACLS  } from '../../../modules/nf-core/custom/tabulartogseacls/main.nf'
include { CUSTOM_TABULARTOGSEACHIP } from '../../../modules/nf-core/custom/tabulartogseachip/main.nf'
include { GSEA_GSEA                } from '../../../modules/nf-core/gsea/gsea/main.nf'
include { PROPR_GREA               } from "../../../modules/nf-core/propr/grea/main.nf"

// Combine meta maps, including merging non-identical values of shared keys (e.g. 'id')
def mergeMaps(meta, meta2){
    (meta + meta2).collectEntries { k, v ->
        meta[k] && meta[k] != v ? [k, "${meta[k]}_${v}"] : [k, v]
    }
}

workflow DIFFERENTIAL_FUNCTIONAL_ENRICHMENT {
    take:
    // input data for functional analysis
    // Note that genesets and background are optional depending on the method.
    // Please set to [] if not provided, eg: [meta, input, [], [], method]
    ch_input                 // [ meta, input file, genesets file, background file, method to run ]

    // other - for the moment these files are only needed for GSEA
    // as it is the only one that takes expression data as input
    // if in the future this setting is changed, this section could be removed
    ch_contrasts             // [ meta, [meta_contrast], [variable], [reference], [target], [formula], [comparison] ]
    ch_samplesheet           // [ meta, samples sheet ]
    ch_featuresheet          // [ meta, features sheet, features id, features symbol ]

    main:

    ch_versions = Channel.empty()

    // Add method information into meta map of ch_input
    // This information is used later to determine which method to run for each input
    // Also, reorganize the structure to match them with the modules' input organization

    ch_input_for_other = ch_input
        .multiMap {
            meta, file, genesets, background, method ->
            def meta_with_method = meta + [ 'functional_method': method ]
            input:
                [ meta_with_method, file ]
            genesets:
                [ meta_with_method, genesets ]
            background:
                [ meta_with_method, background ]
        }

    // In the case of GSEA, it needs additional files coming from other channels that other methods don't use
    // here we define the input channel for the GSEA section

    def criteria = multiMapCriteria { meta, input, genesets, background, method, samplesheet, featuresheet, features_id, features_symbol, meta_contrast, variable, reference, target, formula, comparison ->
        def meta_with_method = meta + [ 'functional_method': method ]
        def meta_with_contrast = mergeMaps(meta_contrast, meta_with_method)
        input:
            [ meta_with_contrast, input ]
        input_without_contrast:
            [ meta_with_method, input ]
        genesets:
            [ meta_with_contrast, genesets ]
        contrasts_and_samples:
            [ meta_with_contrast, samplesheet ]
        features:
            [ [meta_with_method, featuresheet], [features_id, features_symbol] ]
        meta_key:
            [ meta_with_method, meta_with_contrast ]
    }

    // GSEA uses meta.variable, so only keep contrasts where meta.variable is present
    ch_contrasts_transposed = ch_contrasts.transpose()
        .filter { meta, contrastMap, variable, reference, target, formula, comparison ->
            variable?.trim()
        }

    ch_input_for_gsea = ch_input
        .filter{ it[4] == 'gsea' }
        .combine(ch_samplesheet.join(ch_featuresheet), by:0)
        .combine(ch_contrasts_transposed, by:0)
        .multiMap(criteria)

    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    GPROFILER2_GOST(
        ch_input_for_other.input.filter{ it[0].functional_method == 'gprofiler2' },
        ch_input_for_other.genesets.filter{ it[0].functional_method == 'gprofiler2'},
        ch_input_for_other.background.filter{ it[0].functional_method == 'gprofiler2'}
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    // CLS input can be as many as combinations of input x contrasts
    // But not for GCT and CHIP

    CUSTOM_TABULARTOGSEACLS( ch_input_for_gsea.contrasts_and_samples )

    CUSTOM_TABULARTOGSEACHIP(
        ch_input_for_gsea.features.unique().map{it[0]},
        ch_input_for_gsea.features.unique().map{it[1]}
    )

    CUSTOM_TABULARTOGSEAGCT( ch_input_for_gsea.input_without_contrast.unique() )

    // switch the meta by the one with contrast
    ch_chip_gct = CUSTOM_TABULARTOGSEACHIP.out.chip
        .join( CUSTOM_TABULARTOGSEAGCT.out.gct )
        .combine( ch_input_for_gsea.meta_key, by:0 )
        .map { meta, chip, gct, meta_with_contrast ->
            return [ meta_with_contrast, chip, gct ]
        }

    // prepare the input for GSEA
    ch_input_for_gsea = ch_chip_gct
        .join( CUSTOM_TABULARTOGSEACLS.out.cls )
        .join( ch_input_for_gsea.genesets )
        .multiMap { meta, chip, gct, cls, genesets ->
            input:
                [ meta, gct, cls, genesets ]
            contrast:
                [ meta.reference, meta.target ]
            chip:
                [ meta, chip ]
        }

    GSEA_GSEA(
        ch_input_for_gsea.input,
        ch_input_for_gsea.contrast,
        ch_input_for_gsea.chip
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    PROPR_GREA(
        ch_input_for_other.input.filter{ it[0].functional_method == 'grea' },
        ch_input_for_other.genesets.filter{ it[0].functional_method == 'grea' }
    )

    // collect versions info
    ch_versions = ch_versions
        .mix(GPROFILER2_GOST.out.versions)
        .mix(CUSTOM_TABULARTOGSEAGCT.out.versions)
        .mix(CUSTOM_TABULARTOGSEACLS.out.versions)
        .mix(CUSTOM_TABULARTOGSEACHIP.out.versions)
        .mix(GSEA_GSEA.out.versions)
        .mix(PROPR_GREA.out.versions)

    emit:
    // here we emit the outputs that will be useful afterwards in the
    // nf-core/differentialabundance pipeline

    // gprofiler2-specific outputs
    gprofiler2_all_enrich = GPROFILER2_GOST.out.all_enrich
    gprofiler2_sub_enrich = GPROFILER2_GOST.out.sub_enrich
    gprofiler2_plot_html  = GPROFILER2_GOST.out.plot_html

    // gsea-specific outputs
    gsea_report           = GSEA_GSEA.out.report_tsvs_ref.join(GSEA_GSEA.out.report_tsvs_target)

    // grea-specific outputs
    grea_results          = PROPR_GREA.out.results

    // tool versions
    versions              = ch_versions
}
