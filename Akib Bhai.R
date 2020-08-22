akib<-read_xlsx("G:/ML0151/HSE.xlsx",col_names = T)
# What does the data look like?


# Changing column names to shorter names
colnames(HKkids) <- c("Index", "Height", "Weight")

# Let's check general information  about the data!
str(akib)
akib2<-akib[,c(4,30)]
head(akib2)


# It is always good practice to visualize your data!

# install the ggplot2 package
library(ggplot2)

ggplot(akib, aes(x = akib$H04, y = akib$E02)) +
  geom_point(shape=1, size = 5, color = "black", alpha = 1/3) +
  geom_point(size = 5, color = "red4", alpha = 1/3) +
  labs (x = "Hearing Power", y = "Sound Intensity")

# Calculate and store the mean and standard deviation for each dataset, before we normalize

# "HKkids" is our data
# "centers" is the number of clusters
# "iter.max" is the number of iterations
# Try experimenting with different numbers!
clusters_numbers = kmeans (akib, centers = 3, iter.max = 10)
# We are using a color_offset = 4 just for making things easier to see.
# Try experimenting with different numbers!
color_offset = 4

# We have to convert numbers to categorical data in order to color the chart
akib$Cluster <- as.factor(clusters_numbers$cluster + color_offset)


ggplot() +
  geom_point(data = akib, aes(x = akib$H04, y = akib$E02), shape=5, size = 5, color = "black", alpha = 1/2) +
  geom_point(data = akib, aes(x = akib$H04, y = akib$E02), shape=5, size = 5, colour = akib$Cluster, alpha = 1/2) +
  labs (x = "Hearing Power", y = "Sound Intensity")

