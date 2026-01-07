## Recommended Validation Experiment

### Objective
Test if computational breathing dynamics scores correlate with 
actual CRISPR off-target rates in human cells.

### Cell Line
- HEK293T cells (ATCC CRL-3216)
- Stably expressing SpCas9 (lentiCRISPRv2)

### Guide Selection (n=50)
- 20 "high risk" guides (predicted d>2.0, seed GC changes)
- 20 "low risk" guides (predicted d<1.0, PAM region changes)  
- 10 random controls

### Experimental Protocol
1. Synthesize guide RNAs (order from IDT or Synthego)
2. Transfect cells with gRNA + Cas9
3. Measure outcomes at 48h:
   - On-target editing: T7EE1 assay or NGS
   - Off-targets: GUIDE-seq or targeted deep sequencing
4. Correlate breathing scores with measured off-target rates

### Expected Result
Spearman correlation r>0.6 between breathing score 
and actual off-target frequency validates model.

### Estimated Cost
- Guide synthesis: $50/guide × 50 = $2,500
- Cell culture + reagents: $5,000
- Sequencing: $15,000
- Total: ~$25,000

### Timeline
3-4 months
