library(syuzhet)
my_example_text <- "I begin this story with a neutral statement.  
  Basically this is a very silly test.  
  You are testing the Syuzhet package using short, inane sentences.  
  I am actually very happy today. 
  I have finally finished writing this package.  
  Tomorrow I will be very sad. 
  I won't have anything left to do. 
  I might get angry and decide to do something horrible.  
  I might destroy the entire package and start from scratch.  
  Then again, I might find it satisfying to have completed my first R package. 
  Honestly this use of the Fourier transformation is really quite elegant.  
  You might even say it's beautiful!"
s_v <- get_sentences(my_example_text)

class(s_v)
str(s_v)

poa_word_v <- get_tokens(joyces_portrait, pattern = "\\W")
syuzhet_vector <- get_sentiment(s_v, method="syuzhet")
bing_vector <- get_sentiment(poa_v, method="bing")
afinn_vector <- get_sentiment(poa_v, method="afinn")
nrc_vector <- get_sentiment(poa_v, method="nrc", lang = "english")

rbind(
  sign(head(syuzhet_vector)),
  sign(head(bing_vector)),
  sign(head(afinn_vector)),
  sign(head(nrc_vector))
)

sum(syuzhet_vector)

mean(syuzhet_vector)

summary(syuzhet_vector)

my_example_text <- "I begin this story with a neutral statement.  
  Basically this is a very silly test.  
  You are testing the Syuzhet package using short, inane sentences.  
  I am actually very happy today. 
  I have finally finished writing this package.  
  Tomorrow I will be very sad. 
  I won't have anything left to do. 
  I might get angry and decide to do something horrible.  
  I might destroy the entire package and start from scratch.  
  Then again, I might find it satisfying to have completed my first R package. 
  Honestly this use of the Fourier transformation is really quite elegant.  
  You might even say it's beautiful!"
s_v <- get_sentences(my_example_text)
s_v_sentiment <- get_sentiment(s_v)
plot(
  s_v_sentiment, 
  type="l", 
  main="Example Plot Trajectory", 
  xlab = "Narrative Time", 
  ylab= "Emotional Valence"
  )

plot(
  syuzhet_vector, 
  type="h", 
  main="Example Plot Trajectory", 
  xlab = "Narrative Time", 
  ylab= "Emotional Valence"
)

percent_vals <- get_percentage_values(syuzhet_vector, bins = 10)
plot(
  percent_vals, 
  type="l", 
  main="Joyce's Portrait Using Percentage-Based Means", 
  xlab = "Narrative Time", 
  ylab= "Emotional Valence", 
  col="red"
  )

percent_vals <- get_percentage_values(syuzhet_vector, bins = 20)
plot(
  percent_vals, 
  type="l", 
  main="Joyce's Portrait Using Percentage-Based Means", 
  xlab = "Narrative Time", 
  ylab= "Emotional Valence", 
  col="red"
  )

library(syuzhet)
ft_values <- get_transformed_values(
      syuzhet_vector, 
      low_pass_size = 3, 
      x_reverse_len = 100,
      padding_factor = 2,
      scale_vals = TRUE,
      scale_range = FALSE
      )
plot(
  ft_values, 
  type ="l", 
  main ="Joyce's Portrait using Transformed Values", 
  xlab = "Narrative Time", 
  ylab = "Emotional Valence", 
  col = "red"
  )

library(syuzhet)
dct_values <- get_dct_transform(
      syuzhet_vector, 
      low_pass_size = 5, 
      x_reverse_len = 100,
      scale_vals = F,
      scale_range = T
      )
plot(
  dct_values, 
  type ="l", 
  main ="Joyce's Portrait using Transformed Values", 
  xlab = "Narrative Time", 
  ylab = "Emotional Valence", 
  col = "red"
  )

path_to_a_text_file <- system.file("extdata", "bovary.txt", package = "syuzhet")
bovary <- get_text_as_string(path_to_a_text_file)
bovary_v <- get_sentences(bovary)
bovary_sentiment <- get_sentiment(bovary_v)
simple_plot(bovary_sentiment)

nrc_data <- get_nrc_sentiment(s_v)

angry_items <- which(nrc_data$anger > 0)
s_v[angry_items]

joy_items <- which(nrc_data$joy > 0)
s_v[joy_items]


fear_items <- which(nrc_data$fear > 0)
s_v[fear_items]

install.packages("pander")

library(pander)

pander::pandoc.table(nrc_data[, 1:8], split.table = Inf)
pander::pandoc.table(nrc_data[, 9:10])

valence <- (nrc_data[, 9]*-1) + nrc_data[, 10]
valence
barplot(
  sort(colSums(prop.table(nrc_data[, 1:8]))), 
  horiz = TRUE, 
  cex.names = 0.7, 
  las = 1, 
  main = "Emotions in Sample text", xlab="Percentage"
  )

path_to_a_text_file <- system.file("extdata", "portrait.txt",package = "syuzhet")
joyces_portrait <- get_text_as_string(path_to_a_text_file)
poa_v <- get_sentences(joyces_portrait)
poa_values <- get_sentiment(poa_v, method="syuzhet")

path_to_a_text_file <- system.file("extdata", "bovary.txt", package = "syuzhet")
bovary <- get_text_as_string(path_to_a_text_file)
bovary_v <- get_sentences(bovary)
bovary_values <- get_sentiment(bovary_v)


pwdw <- round(length(poa_values)*.1)
poa_rolled <- zoo::rollmean(poa_values, k=pwdw)
bwdw <- round(length(bovary_values)*.1)
bov_rolled <- zoo::rollmean(bovary_values, k=bwdw)

poa_list <- rescale_x_2(poa_rolled)
bov_list <- rescale_x_2(bov_rolled)

plot(poa_list$x, 
     poa_list$z, 
     type="l", 
     col="blue", 
     xlab="Narrative Time", 
     ylab="Emotional Valence")
lines(bov_list$x, bov_list$z, col="red")

poa_sample <- seq(1, length(poa_list$x), by=round(length(poa_list$x)/100))
bov_sample <- seq(1, length(bov_list$x), by=round(length(bov_list$x)/100))

plot(poa_list$x[poa_sample], 
     poa_list$z[poa_sample], 
     type="l", 
     col="blue",
     xlab="Narrative Time (sampled)", 
     ylab="Emotional Valence"
     )
lines(bov_list$x[bov_sample], bov_list$z[bov_sample], col="red")

dist(rbind(poa_list$z[poa_sample], bov_list$z[bov_sample]))
cor(cbind(poa_list$z[poa_sample], bov_list$z[bov_sample]))

poa_x <- 1:length(poa_values)
poa_y <- poa_values
raw_poa <- loess(poa_y ~ poa_x, span=.5)
poa_line <- rescale(predict(raw_poa))
bov_x <- 1:length(bovary_values)
bov_y <- bovary_values
raw_bov <- loess(bov_y ~ bov_x, span=.5)
bov_line <- rescale(predict(raw_bov))
poa_sample <- seq(1, length(poa_line), by=round(length(poa_line)/100))
bov_sample <- seq(1, length(bov_line), by=round(length(bov_line)/100))
plot(poa_line[poa_sample], 
     type="l", 
     col="blue",
     xlab="Narrative Time (sampled)", 
     ylab="Emotional Valence"
     )
lines(bov_line[bov_sample], col="red")


path_to_a_text_file <- system.file("extdata", "quijote.txt",package = "syuzhet")
my_text <- get_text_as_string(path_to_a_text_file)
char_v <- get_sentences(my_text)
method <- "nrc"
lang <- "spanish"
my_text_values <- get_sentiment(char_v, method=method, language=lang)
my_text_values[1:10]
rm(my_text_values)












































