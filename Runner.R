source('[Analysis] GBM Testing Library.R')

## NOTE: file patient_SP19_platypus was renamed to patient_A23_platypus because that's thes same patient
## NOTE: file A23_SP19.platypus was renamed to A23.platypus because that's thes same patient
## NOTE: file SP28_R11.platypus was renamed to SP28.platypus because that's thes same patient
patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A34', 'A44', 'SP28')


setwd('[Data] CCFs/')
source('[Analysis] 1.Extractor.R')
source('[Analysis] 2.Extract Clonal SNVs.R')
source('[Analysis] 3.Extract Exome SNVs.R')


setwd('../[Data] TES_1')
source('[Analysis] 1.Testable.R')
source('[Analysis] 2.Margin-S-readcountss.R')

setwd('../[Data] TES_2')
source('[Analysis] 1.Testable.R')
source('[Analysis] 2.Margin-S-readcountss.R')


setwd('../[Data] WES_PASS')
source('[Analysis] 1.Extractor.R')

# A34 has not SNV to test
patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A44', 'SP28')

source('[Analysis] 2.Beta-Binomial Fit.R')

patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A34', 'A44', 'SP28')

source('[Analysis] 3.VAFs.R')
source('[Analysis] 4.Exome-SNVs.R')
source('[Analysis] 5.Margin-S-readcountss.R')
# source('[Analysis] 6.Gene-Mut-Annotations.R')


setwd('..')
source('[Analysis] GBM Testing H0.R')

