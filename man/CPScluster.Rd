\name{CPScluster}
\alias{CPScluster}
\docType{data}
\title{
CPS subset with clusterings
}
\description{
Different hierarchical clusterings and k-means clusterings as well as a model-based clustering have been applied to several financial variables for a random sample of ten thousand observations.
}
\usage{data(CPScluster)}
\format{
  A data frame with 10000 observations on the following 39 variables.
  \describe{
    \item{\code{Age}}{a numeric vector}
    \item{\code{Sex}}{a factor with levels \code{female} \code{male}}
    \item{\code{Race}}{a factor with levels \code{Black} \code{White}}
    \item{\code{Ethnic}}{a factor with levels \code{Central and South American} \code{Chicano} \code{Cuban} \code{Mexican American} \code{Mexicano} \code{Other Spanish} \code{Others} \code{Puerto Rican} \code{don t know}}
    \item{\code{Marital.Status}}{a factor with levels \code{Divorced} \code{Married} \code{Never married} \code{Separated} \code{Widowed}}
    \item{\code{Kind.of.Family}}{a factor with levels \code{Husband-Wife family} \code{Other female head family} \code{Other male head family}}
    \item{\code{Classical}}{a factor with levels \code{All other} \code{Classical Husband-Wife family}}
    \item{\code{Family.Type}}{a factor with levels \code{Nonfamily householder} \code{Primary family} \code{Related subfamily} \code{Unrelated individual} \code{Unrelated subfamily}}
    \item{\code{Number.of.Persons.in.Family}}{a numeric vector}
    \item{\code{Number.of.Kids}}{a numeric vector}
    \item{\code{Education.of.Head}}{a factor with levels \code{00. Children} \code{01. Less than 1st grade} \code{02. 1st, 2nd, 3rd or 4th grade} \code{03. 5th or 6th grade} \code{04. 7th or 8th grade} \code{05. 9th grade} \code{06. 10th grade} \code{07. 11th grade} \code{08. 12th grade, no diploma} \code{09. High school graduate, high school diploma or equivalent} \code{10. Some college but no degree} \code{11. Ass. college deg.:occup./voc. program} \code{12. Ass. college deg.:academic program} \code{13. Bachelors degree} \code{14. Masters degree} \code{15. Professional school degree} \code{16. Doctoral degree}}
    \item{\code{Labor.Status}}{a factor with levels \code{Armed forces} \code{Children} \code{Disabled} \code{Employed} \code{Employed, not at work} \code{Not in labor force, other} \code{Retired} \code{Unemployed, looking for work} \code{Unemployed, on layoff}}
    \item{\code{Class.of.Worker}}{a factor with levels \code{Federal government} \code{Local government} \code{Not in labor force} \code{Private sector} \code{Self-empl.: incorporated} \code{Self-empl.: unincorporated} \code{State government} \code{Unemployed, no previous experience} \code{Without pay}}
    \item{\code{Working.Hours}}{a numeric vector}
    \item{\code{Income.of.Head}}{a numeric vector}
    \item{\code{Family.Income}}{a numeric vector}
    \item{\code{Taxable.Income}}{a numeric vector}
    \item{\code{Federal.tax}}{a numeric vector}
    \item{\code{Family.sequence.number}}{a numeric vector}
    \item{\code{State}}{a factor with levels \code{Alabama} \code{Alaska} \code{Arizona} \code{Arkansas} \code{California} \code{Colorado} \code{Connecticut} \code{Delaware} \code{District of Columbia} \code{Florida} \code{Georgia} \code{Hawaii} \code{Idaho} \code{Illinois} \code{Indiana} \code{Iowa} \code{Kansas} \code{Kentucky} \code{Louisiana} \code{Maine} \code{Maryland} \code{Massachusetts} \code{Michigan} \code{Minnesota} \code{Mississippi} \code{Missouri} \code{Montana} \code{Nebraska} \code{Nevada} \code{New Hampshire} \code{New Jersey} \code{New Mexico} \code{New York} \code{North Carolina} \code{North Dakota} \code{Ohio} \code{Oklahoma} \code{Oregon} \code{Pennsylvania} \code{Rhode Island} \code{South Carolina} \code{South Dakota} \code{Tennessee} \code{Texas} \code{Utah} \code{Vermont} \code{Virginia} \code{Washington} \code{West Virginia} \code{Wisconsin} \code{Wyoming}}
    \item{\code{Division}}{a factor with levels \code{East North Central} \code{East South Central} \code{Middle Atlantic} \code{Mountain} \code{New England } \code{Pacific} \code{South Atlantic} \code{West North Central} \code{West South Central}}
    \item{\code{Region}}{a factor with levels \code{Midwest} \code{North East} \code{South} \code{West}}
    \item{\code{hc4}}{a numeric vector}
    \item{\code{hc6}}{a numeric vector}
    \item{\code{hc8}}{a numeric vector}
    \item{\code{hc12}}{a numeric vector}
    \item{\code{hcs4}}{a numeric vector}
    \item{\code{hcs6}}{a numeric vector}
    \item{\code{hcs8}}{a numeric vector}
    \item{\code{hcs12}}{a numeric vector}
    \item{\code{hcw4}}{a numeric vector}
    \item{\code{hcw6}}{a numeric vector}
    \item{\code{hcw8}}{a numeric vector}
    \item{\code{hcw12}}{a numeric vector}
    \item{\code{km4}}{a numeric vector}
    \item{\code{km6}}{a numeric vector}
    \item{\code{km8}}{a numeric vector}
    \item{\code{km12}}{a numeric vector}
    \item{\code{mc12}}{a numeric vector}
  }
}

\keyword{datasets}
