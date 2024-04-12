install.packages('mediation')
library(mediation)
packageVersion('mediation')


data <- read.csv('/proj/gibbons/kramos/github/IBT-and-the-Gut-Microbiome/american_gut_cohort/american_gut_filtered_standardized.csv')

# Make sex a factor variable
data$sex <- factor(data$sex)

# Create my mediator model 
model.MH <- lm(formula=simpson~age+sex+bowel_movement_frequency+BMI_CALC+high_vegetable_consumption+height_cm, data=data)
summary(model.MH)

# Create model.Y
model.YH <- glm(formula=has_cdiff~age+sex+bowel_movement_frequency+BMI_CALC+height_cm+high_vegetable_consumption+simpson, data=data)

resultsH <- mediate(model.MH, model.YH, treat='height_cm', mediator='simpson', boot=TRUE, sims=5001)
summary(resultsH)

################ Do this again for diet ################

#Create my mediator model 
model.MV <- lm(formula=simpson~age+sex+bowel_movement_frequency+BMI_CALC+high_vegetable_consumption, data=data)
summary(model.MV)

#Create model.Y
model.YV <- glm(formula=has_cdiff~age+sex+bowel_movement_frequency+BMI_CALC+high_vegetable_consumption+simpson, data=data)

resultsV<-mediate(model.MV, model.YV, treat='high_vegetable_consumption', mediator='simpson', boot=TRUE, sims=5001)
summary(resultsV)
