### Importing the data into dat: ###
dat = read.csv("W:/somnath.datta/shoumisarkar/Fluorosis/RealData/Dental Exam Data/permanent_fluorosis_hypoplasia_opacity_zone.txt",
           sep="", stringsAsFactors=TRUE)


library(tidyverse)

#subset by ages:
age23 = dat %>% select(SUBJECT_ID, Zone, FRI_Score23)
age17 = dat %>% select(SUBJECT_ID, Zone, FRI_Score17)
age13 = dat %>% select(SUBJECT_ID, Zone, FRI_Score13)
age9 = dat %>% select(SUBJECT_ID, Zone, FRI_Score9)

### Age 9 ###

plot9 = 
  ggplot(age9) + 
  aes(FRI_Score9, fill=Zone) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), position="dodge")  + #side by side
  labs(title = "Age 9", x = "FRI", y = "Proportions")

### Age 13 ###

plot13 =
  ggplot(age13) + 
  aes(FRI_Score13, fill=Zone) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), position="dodge")  + #side by side
  labs(title = "Age 13", x = "FRI", y = "Proportions")

### Age 17 ###

plot17 =
  ggplot(age17) + 
  aes(FRI_Score17, fill=Zone) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), position="dodge")  + #side by side
  labs(title = "Age 17", x = "FRI", y = "Proportions")


### Age 23 ###

plot23 =
  ggplot(age23) + 
  aes(FRI_Score23, fill=Zone) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), position="dodge")  + #side by side
  labs(title = "Age 23", x = "FRI", y = "Proportions")


require(gridExtra)
grid.arrange(plot9, plot13, plot17, plot23, ncol=4)


####################################################################################
############### Convert to Fluorosis Present/Absent ################################
####################################################################################

#subset by ages:
age23 = dat %>% select(SUBJECT_ID, Zone, FRI_Score23)
age17 = dat %>% select(SUBJECT_ID, Zone, FRI_Score17)
age13 = dat %>% select(SUBJECT_ID, Zone, FRI_Score13)
age9 = dat %>% select(SUBJECT_ID, Zone, FRI_Score9)

age9$FRI_Score9 = ifelse(age9$FRI_Score9==0, 0,1)
age13$FRI_Score13 = ifelse(age13$FRI_Score13==0, 0,1)
age17$FRI_Score17 = ifelse(age17$FRI_Score17==0, 0,1)
age23$FRI_Score23 = ifelse(age23$FRI_Score23==0, 0,1)


### Age 9 ###

plot9 <- ggplot(age9) + 
  aes(FRI_Score9, fill = Zone) + 
  geom_bar(aes(y = (..count..) / sum(..count..)), position = "dodge") +
  labs(title = "Age 9", x = "FRI", y = "Proportions") +
  scale_x_continuous(breaks = c(0, 1))  # Set the breaks to 0 and 1 on the x-axis

# Print the plot
print(plot9)

### Age 13 ###

plot13 =
  ggplot(age13) + 
  aes(FRI_Score13, fill=Zone) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), position="dodge")  + #side by side
  labs(title = "Age 13", x = "FRI", y = "Proportions") +
  scale_x_continuous(breaks = c(0, 1))  # Set the breaks to 0 and 1 on the x-axis

### Age 17 ###

plot17 =
  ggplot(age17) + 
  aes(FRI_Score17, fill=Zone) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), position="dodge")  + #side by side
  labs(title = "Age 17", x = "FRI", y = "Proportions") +
  scale_x_continuous(breaks = c(0, 1))  # Set the breaks to 0 and 1 on the x-axis


### Age 23 ###

plot23 =
  ggplot(age23) + 
  aes(FRI_Score23, fill=Zone) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), position="dodge")  + #side by side
  labs(title = "Age 23", x = "FRI", y = "Proportions") +
  scale_x_continuous(breaks = c(0, 1))  # Set the breaks to 0 and 1 on the x-axis


require(gridExtra)
grid.arrange(plot9, plot13, plot17, plot23, ncol=4)

