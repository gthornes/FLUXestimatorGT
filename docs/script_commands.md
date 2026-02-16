# Example command line inputs for each script

## scFEA script:

```
$env:KMP_DUPLICATE_LIB_OK="TRUE"; cd Z:\thornes\Pavlicev_lab_rotation\scFEA; python src\scFEA.py --data_dir data --input_dir ../FLUXestimator/results/scfea/combined.29.01.26 --res_dir ../FLUXestimator/results/scfea/combined.29.01.26 --test_file expression_input.csv --moduleGene_file module_gene_complete_mouse_m168.csv
```

## Visualisation script:

```
python scripts\visualisation.py -i results\scfea\22.01.26\expression_input_module168_cell7964_batch7964_LR0.008_epoch100_SCimpute_F_lambBal1_lambSca1_lambCellCor1_lambModCor_1e-2_20260122-160247.csv -o results\figures\22.01.26
```

## Annotate modules script:

```
python scripts\annotate_modules.py --flux results\scfea\combined.26.01.26\expression_input_module168_cell10779_batch10779_LR0.008_epoch100_SCimpute_F_lambBal1_lambSca1_lambCellCor1_lambModCor_1e-2_20260127-100331.csv -o results\tables\combined.26.01.26 --create-annotated-file
```

## Statistical analysis script

```
python scripts\statistical_analysis.py -i results\scfea\combined.29.01.26\expression_input_module168_cell9160_batch9160_LR0.008_epoch100_SCimpute_F_lambBal1_lambSca1_lambCellCor1_lambModCor_1e-2_20260130-111417_annotated.csv -o results\statistical_analysis\combined.29.01.26 --fdr 0.15
```