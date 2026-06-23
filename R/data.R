#' FRED-MD Macroeconomic Dataset for Inflation Forecasting
#'
#' Monthly macroeconomic and financial variables from the FRED-MD database
#' (McCracken and Ng, 2016), prepared for inflation forecasting.
#' Variables are transformed following the FRED-MD recommended transformations
#' (stationarity-inducing), except UNRATE which is kept in levels as a
#' near-unit-root process. The target variable \code{infl_cpi} is constructed
#' as \eqn{100 \times \Delta \log(\text{CPIAUCSL})}{100 * diff(log(CPIAUCSL))}.
#' CPIAUCSL itself is excluded from the predictor set.
#'
#' @format A tibble with 785 rows and 113 variables:
#' \describe{
#'   \item{date}{observation date (monthly, 1959-12 to 2025-04)}
#'   \item{infl_cpi}{CPI inflation rate (target variable)}
#'   \item{UNRATE}{unemployment rate (levels, not transformed)}
#'   \item{RPI}{real personal income}
#'   \item{W875RX1}{real personal income ex transfer receipts}
#'   \item{DPCERA3M086SBEA}{real personal consumption expenditures}
#'   \item{CMRMTSPLx}{real mfg and trade industries sales}
#'   \item{RETAILx}{retail and food services sales}
#'   \item{INDPRO}{IP index}
#'   \item{IPFPNSS}{IP final products and nonindustrial supplies}
#'   \item{IPFINAL}{IP final products}
#'   \item{IPCONGD}{IP consumer goods}
#'   \item{IPDCONGD}{IP durable consumer goods}
#'   \item{IPNCONGD}{IP nondurable consumer goods}
#'   \item{IPBUSEQ}{IP business equipment}
#'   \item{IPMAT}{IP materials}
#'   \item{IPDMAT}{IP durable materials}
#'   \item{IPNMAT}{IP nondurable materials}
#'   \item{IPMANSICS}{IP manufacturing SIC}
#'   \item{IPFUELS}{IP fuels}
#'   \item{CUMFNS}{capacity utilization}
#'   \item{HWI}{help-wanted index}
#'   \item{HWIURATIO}{ratio help-wanted/unemployment}
#'   \item{CLF16OV}{civilian labor force}
#'   \item{CE16OV}{civilian employment}
#'   \item{UEMPMEAN}{avg duration of unemployment}
#'   \item{UEMPLT5}{unemployment less than 5 weeks}
#'   \item{UEMP5TO14}{unemployment 5-14 weeks}
#'   \item{UEMP15OV}{unemployment 15+ weeks}
#'   \item{UEMP15T26}{unemployment 15-26 weeks}
#'   \item{UEMP27OV}{unemployment 27+ weeks}
#'   \item{CLAIMSx}{initial claims}
#'   \item{PAYEMS}{all employees total nonfarm}
#'   \item{USGOOD}{all employees goods-producing}
#'   \item{CES1021000001}{all employees mining and logging}
#'   \item{USCONS}{all employees construction}
#'   \item{MANEMP}{all employees manufacturing}
#'   \item{DMANEMP}{all employees durable goods}
#'   \item{NDMANEMP}{all employees nondurable goods}
#'   \item{SRVPRD}{all employees service-providing}
#'   \item{USTPU}{all employees trade transportation utilities}
#'   \item{USWTRADE}{all employees wholesale trade}
#'   \item{USTRADE}{all employees retail trade}
#'   \item{USFIRE}{all employees financial activities}
#'   \item{USGOVT}{all employees government}
#'   \item{CES0600000007}{avg weekly hours goods-producing}
#'   \item{AWOTMAN}{avg weekly overtime hours manufacturing}
#'   \item{AWHMAN}{avg weekly hours manufacturing}
#'   \item{HOUST}{housing starts total}
#'   \item{HOUSTNE}{housing starts northeast}
#'   \item{HOUSTMW}{housing starts midwest}
#'   \item{HOUSTS}{housing starts south}
#'   \item{HOUSTW}{housing starts west}
#'   \item{AMDMNOx}{new orders durable goods}
#'   \item{AMDMUOx}{unfilled orders durable goods}
#'   \item{BUSINVx}{total business inventories}
#'   \item{ISRATIOx}{total business inventories to sales ratio}
#'   \item{M1SL}{M1 money stock}
#'   \item{M2SL}{M2 money stock}
#'   \item{M2REAL}{real M2 money stock}
#'   \item{TOTRESNS}{total reserves}
#'   \item{NONBORRES}{nonborrowed reserves}
#'   \item{BUSLOANS}{commercial and industrial loans}
#'   \item{REALLN}{real estate loans}
#'   \item{NONREVSL}{total nonrevolving credit}
#'   \item{CONSPI}{consumer sentiment index}
#'   \item{S&P 500}{S&P 500 index}
#'   \item{S&P div yield}{S&P 500 dividend yield}
#'   \item{S&P PE ratio}{S&P 500 price-earnings ratio}
#'   \item{FEDFUNDS}{effective federal funds rate}
#'   \item{TB3MS}{3-month treasury bill}
#'   \item{TB6MS}{6-month treasury bill}
#'   \item{GS1}{1-year treasury rate}
#'   \item{GS5}{5-year treasury rate}
#'   \item{GS10}{10-year treasury rate}
#'   \item{AAA}{Moody's AAA corporate bond yield}
#'   \item{BAA}{Moody's BAA corporate bond yield}
#'   \item{TB3SMFFM}{3m treasury spread}
#'   \item{TB6SMFFM}{6m treasury spread}
#'   \item{T1YFFM}{1yr treasury spread}
#'   \item{T5YFFM}{5yr treasury spread}
#'   \item{T10YFFM}{10yr treasury spread}
#'   \item{AAAFFM}{AAA-FF spread}
#'   \item{BAAFFM}{BAA-FF spread}
#'   \item{EXSZUSx}{Switzerland/US exchange rate}
#'   \item{EXJPUSx}{Japan/US exchange rate}
#'   \item{EXUSUKx}{US/UK exchange rate}
#'   \item{EXCAUSx}{Canada/US exchange rate}
#'   \item{WPSFD49207}{PPI finished goods}
#'   \item{WPSFD49502}{PPI finished consumer goods}
#'   \item{WPSID61}{PPI intermediate materials}
#'   \item{WPSID62}{PPI crude materials}
#'   \item{OILPRICEx}{crude oil prices}
#'   \item{PPICMM}{PPI metals and metal products}
#'   \item{CPIAPPSL}{CPI apparel}
#'   \item{CPITRNSL}{CPI transportation}
#'   \item{CPIMEDSL}{CPI medical care}
#'   \item{CUSR0000SAC}{CPI commodities}
#'   \item{CUSR0000SAD}{CPI durables}
#'   \item{CUSR0000SAS}{CPI services}
#'   \item{CPIULFSL}{CPI all items less food}
#'   \item{CUSR0000SA0L2}{CPI all items less shelter}
#'   \item{CUSR0000SA0L5}{CPI all items less medical care}
#'   \item{PCEPI}{PCE price index}
#'   \item{DDURRG3M086SBEA}{PCE durable goods}
#'   \item{DNDGRG3M086SBEA}{PCE nondurable goods}
#'   \item{DSERRG3M086SBEA}{PCE services}
#'   \item{CES0600000008}{avg hourly earnings goods-producing}
#'   \item{CES2000000008}{avg hourly earnings construction}
#'   \item{CES3000000008}{avg hourly earnings manufacturing}
#'   \item{DTCOLNVHFNM}{consumer motor vehicle loans}
#'   \item{DTCTHFNM}{total consumer loans}
#'   \item{INVEST}{investment variable}
#' }
#' @source FRED-MD \url{https://www.stlouisfed.org/research/economists/mccracken/fred-databases}
"fredmd"

#' Raw FRED-MD Data (Before Transformation)
#'
#' Untransformed monthly macroeconomic series from the FRED-MD database,
#' with \code{infl_cpi} constructed as \eqn{100 \times \Delta \log(\text{CPIAUCSL})}.
#' CPIAUCSL itself is excluded (same as \code{\link{fredmd}}).
#'
#' @format A tibble with the same column names as \code{\link{fredmd}}
#'   but with untransformed values.
#' @source FRED-MD \url{https://www.stlouisfed.org/research/economists/mccracken/fred-databases}
"data_raw"

#' Predictor Variable Names from FRED-MD
#'
#' Character vector of 111 predictor variable names from the \code{\link{fredmd}}
#' dataset. Includes all columns except \code{date} and the target \code{infl_cpi}.
#' UNRATE is included as a predictor (kept in levels).
#'
#' @format A character vector of length 111.
"vars_all"
