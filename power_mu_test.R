power_mu_tests <- function(sample_size = 10, to = 2, iter = 10^4, signif_level = 0.05){

  # パッケージの読み込み
  library(dplyr)
  library(ggplot2)
  library(gridExtra)

  # mu: 位置母数とする。（注意 : 正規分布の場合は母平均のこと。）
  # 検定問題 H0 : mu == 0 v.s. H1 : mu != 0に対する検定の比較
  # 正規分布の場合と自由度1のt分布（コーシー分布）の場合で、第1種の誤りと第2種の誤りを計算する。
  # ・正規分布 : 一般的な分布
  # ・コーシー分布 : 外れ値が得られやすい
  # 比較する検定
  # (1) 1標本t検定 : t.test関数
  # (2) Wilcoxonの符号付き順位検定 : wilcox.test関数
  # (3) 符号検定 : binom.test関数

  # 引数の説明
  # sample_size : 検出力を調べたいサンプルサイズ
  # to : H1が正しい場合、その具体的なmuの値を0からtoまで調べる。
  # iter : 繰り返し数
  # signif_level : 検定の有意水準

  # 正規分布と自由度1のt分布からサンプリングする。
  mu <- seq(from = from, to = to, by = (to-from)/100)

  sample_norm <- rnorm(101 * sample_size * iter, mean = 0, sd = 1)
  sample_norm <- rep(mu, each = sample_size * iter) + sample_norm
  sample_norm <- matrix(sample_norm, nrow = 101 * iter, byrow = TRUE)

  sample_Cauchy <- rt(101 * sample_size * iter, df = 1)
  sample_Cauchy <- rep(mu, each = sample_size * iter) + sample_Cauchy
  sample_Cauchy <- matrix(sample_Cauchy, nrow = 101 * iter, byrow = TRUE)

  # 検定を実行する。
  res_norm_T <- apply(sample_norm, 1, function(x){t.test(x,mu=0)$p.value <= signif_level})
  res_norm_W <- apply(sample_norm, 1, function(x){wilcox.test(x,mu=0)$p.value <= signif_level})
  res_norm_S <- apply(sample_norm, 1, function(x){binom.test(c(sum(x >= 0), sum(x < 0)))$p.value <= signif_level})
  res_norm <- c(res_norm_T, res_norm_W, res_norm_S)

  res_Cauchy_T <- apply(sample_Cauchy, 1, function(x){t.test(x,mu=0)$p.value <= signif_level})
  res_Cauchy_W <- apply(sample_Cauchy, 1, function(x){wilcox.test(x,mu=0)$p.value <= signif_level})
  res_Cauchy_S <- apply(sample_Cauchy, 1, function(x){binom.test(c(sum(x >= 0), sum(x < 0)))$p.value <= signif_level})
  res_Cauchy <- c(res_Cauchy_T, res_Cauchy_W, res_Cauchy_S)

  res_norm <- apply(matrix(res_norm, ncol = iter, byrow = TRUE), 1, mean)    # 各muの値に対してH0を棄却した確率
  res_Cauchy <- apply(matrix(res_Cauchy, ncol = iter, byrow = TRUE), 1, mean)    # 各muの値に対してH0を棄却した確率

  # 結果をdata_frameにまとめて、ggplot2に流し込む準備をする。
  data <- data_frame(真の位置母数 = rep(rep(mu, times = 3), times = 2),
                           検定 = rep(rep(c("1標本t検定", "Wilcoxonの符号付き順位検定", "符号検定"), each = 101), times = 2),
                           分布 = rep(c("正規分布", "Cauchy分布"), each = 303),
                           棄却率 = c(res_norm, res_Cauchy))

  # 出力を準備する。
  # (1) 標準正規分布とCauchy分布のグラフ
  g1 <- ggplot(data = data_frame(x = c(-4,4)), aes(x = x)) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1))+
    stat_function(fun = dt, args = list(df = 1), color = "red") +
    ylab("確率密度")+
    ggtitle("黒 : 標準正規分布の確率密度関数, 赤 : Cauchy分布の確率密度関数")

  # (2) 正規分布のケース : 第1種の誤りの比較（棒グラフ）
  alpha <- data %>% filter(真の位置母数 == 0 & 分布 == "正規分布")
  g2 <- ggplot(alpha, aes(x = 検定, y = 棄却率))  +
    geom_bar(stat = "identity") +
    xlab("各種検定") +
    ylab("第1種の誤り") +
    ggtitle("各検定の第1種の誤り（正規分布の場合）")

  # (3) 正規分布のケース : 位置母数muに対する各検定の検出力
  power <- data %>% filter(分布 == "正規分布")
  g3 <- ggplot(power, aes(x = 真の位置母数, y = 棄却率, color = 検定)) +
    geom_line() +
    xlab("真の位置母数") +
    ylab("検出力") +
    ggtitle("各検定の検出力（正規分布の場合）")

  # (4) Cauchy分布のケース : 第1種の誤りの比較（棒グラフ）
  alpha <- data %>% filter(mu == 0 & 分布 == "Cauchy分布")
  g4 <- ggplot(alpha, aes(x = 検定, y = 棄却率)) +
    geom_bar(stat = "identity") +
    xlab("各種検定") +
    ylab("第1種の誤り") +
    ggtitle("各検定の第1種の誤り（Cauchy分布の場合）")

  # (5) Cauchy分布のケース : 位置母数muに対する各検定の検出力
  power <- data %>% filter(分布 == "Cauchy分布")
  g5 <- ggplot(power, aes(x = 真の位置母数, y = 棄却率, color = 検定)) +
    geom_line() +
    xlab("真の位置母数") +
    ylab("検出力") +
    ggtitle("各検定の検出力（Cauchy分布の場合）")

  # レイアウトを行列にし、layoutに保存しておく
  layout <- rbind(c(1, 1),
                  c(2, 3),
                  c(4, 5))
  # まとめて1枚に出力
  grid.arrange(g1, g2, g3, g4, g5,
               layout_matrix = layout)

}
