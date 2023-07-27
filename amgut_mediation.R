install.packages('mediation')
library(mediation)

# Import american gut cohort data
data <- read.csv('.../american_gut_cohort/american_gut_filtered_standardized.csv')

# Make sex a factor variable
data$sex <- factor(data$sex)

# Create my mediator model 
model.MH <- lm(formula=simpson~age+sex+bowel_movement_frequency+BMI_CALC+height_cm, data=data)
summary(model.MH)

# Create model.Y
model.YH <- glm(formula=has_cdiff~age+sex+bowel_movement_frequency+BMI_CALC+height_cm+simpson, data=data)

resultsH <- mediate(model.MH, model.YH, treat='height_cm', mediator='simpson', boot=TRUE, sims=5000)
summary(resultsH)


################ Do this again for diet ################

#Create my mediator model 
model.MV <- lm(formula=simpson~age+sex+bowel_movement_frequency+BMI_CALC+high_vegetable_consumption, data=data)
summary(model.MV)

#Create model.Y
model.YV <- glm(formula=has_cdiff~age+sex+bowel_movement_frequency+BMI_CALC+high_vegetable_consumption+simpson, data=data)

resultsV <- mediate(model.MV, model.YV, treat='high_vegetable_consumption', mediator='simpson', boot=TRUE, sims=5000)
summary(resultsV)
