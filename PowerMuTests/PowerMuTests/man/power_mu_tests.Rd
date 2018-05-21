\name{power_mu_tests}
\alias{power_mu_tests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
power_mu_tests(sample_size = 10, from = -2, to = 2, iter = 10^4, signif_level = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_size}{
%%     ~~Describe \code{sample_size} here~~
}
  \item{from}{
%%     ~~Describe \code{from} here~~
}
  \item{to}{
%%     ~~Describe \code{to} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{signif_level}{
%%     ~~Describe \code{signif_level} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (sample_size = 10, from = -2, to = 2, iter = 10^4, signif_level = 0.05) 
{
    library(dplyr)
    library(ggplot2)
    library(gridExtra)
    if (from > 0 | to < 0) {
        print("エラー : from引数は負の数しか、to引数は正の数しか代入できません。")
    }
    mu <- seq(from = from, to = to, by = (to - from)/100)
    sample_norm <- rnorm(101 * sample_size * iter, mean = 0, 
        sd = 1)
    sample_norm <- rep(mu, each = sample_size * iter) + sample_norm
    sample_norm <- matrix(sample_norm, nrow = 101 * iter, byrow = TRUE)
    sample_Cauchy <- rt(101 * sample_size * iter, df = 1)
    sample_Cauchy <- rep(mu, each = sample_size * iter) + sample_Cauchy
    sample_Cauchy <- matrix(sample_Cauchy, nrow = 101 * iter, 
        byrow = TRUE)
    res_norm_T <- apply(sample_norm, 1, function(x) {
        t.test(x, mu = 0)$p.value <= signif_level
    })
    res_norm_W <- apply(sample_norm, 1, function(x) {
        wilcox.test(x, mu = 0)$p.value <= signif_level
    })
    res_norm_S <- apply(sample_norm, 1, function(x) {
        binom.test(c(sum(x >= 0), sum(x < 0)))$p.value <= signif_level
    })
    res_norm <- c(res_norm_T, res_norm_W, res_norm_S)
    res_Cauchy_T <- apply(sample_Cauchy, 1, function(x) {
        t.test(x, mu = 0)$p.value <= signif_level
    })
    res_Cauchy_W <- apply(sample_Cauchy, 1, function(x) {
        wilcox.test(x, mu = 0)$p.value <= signif_level
    })
    res_Cauchy_S <- apply(sample_Cauchy, 1, function(x) {
        binom.test(c(sum(x >= 0), sum(x < 0)))$p.value <= signif_level
    })
    res_Cauchy <- c(res_Cauchy_T, res_Cauchy_W, res_Cauchy_S)
    res_norm <- apply(matrix(res_norm, ncol = iter, byrow = TRUE), 
        1, mean)
    res_Cauchy <- apply(matrix(res_Cauchy, ncol = iter, byrow = TRUE), 
        1, mean)
    data <- data_frame(真の位置母数 = rep(rep(mu, times = 3), 
        times = 2), 検定 = rep(rep(c("1標本t検定", "Wilcoxonの符号付き順位検定", 
        "符号検定"), each = 101), times = 2), 分布 = rep(c("正規分布", 
        "Cauchy分布"), each = 303), 棄却率 = c(res_norm, res_Cauchy))
    g1 <- ggplot(data = data_frame(x = c(-4, 4)), aes(x = x)) + 
        stat_function(fun = dnorm, args = list(mean = 0, sd = 1)) + 
        stat_function(fun = dt, args = list(df = 1), color = "red") + 
        ylab("確率密度") + ggtitle("黒 : 標準正規分布の確率密度関数, 赤 : Cauchy分布の確率密度関数")
    alpha <- data \%>\% filter(真の位置母数 == 0 & 分布 == "正規分布")
    g2 <- ggplot(alpha, aes(x = 検定, y = 棄却率)) + geom_bar(stat = "identity") + 
        xlab("各種検定") + ylab("第1種の誤り") + ggtitle("各検定の第1種の誤り（正規分布の場合）")
    power <- data \%>\% filter(分布 == "正規分布")
    g3 <- ggplot(power, aes(x = 真の位置母数, y = 棄却率, color = 検定)) + 
        geom_line() + xlab("真の位置母数") + ylab("検出力") + 
        ggtitle("各検定の検出力（正規分布の場合）")
    alpha <- data \%>\% filter(mu == 0 & 分布 == "Cauchy分布")
    g4 <- ggplot(alpha, aes(x = 検定, y = 棄却率)) + geom_bar(stat = "identity") + 
        xlab("各種検定") + ylab("第1種の誤り") + ggtitle("各検定の第1種の誤り（Cauchy分布の場合）")
    power <- data \%>\% filter(分布 == "Cauchy分布")
    g5 <- ggplot(power, aes(x = 真の位置母数, y = 棄却率, color = 検定)) + 
        geom_line() + xlab("真の位置母数") + ylab("検出力") + 
        ggtitle("各検定の検出力（Cauchy分布の場合）")
    layout <- rbind(c(1, 1), c(2, 3), c(4, 5))
    grid.arrange(g1, g2, g3, g4, g5, layout_matrix = layout)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
