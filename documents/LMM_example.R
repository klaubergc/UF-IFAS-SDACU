###############################################################
# Workshop: Linear Mixed Models (LMM) in R
# Example: Tree growth of evergreen vs. deciduous species
###############################################################

# Research Question:
# Do evergreen species grow more or less than deciduous species,
# while accounting for tree height and species-level variation?

# Hypotheses:
# H0: Evergreen and deciduous species show no difference in growth.
# H1: Evergreen and deciduous species differ in growth rates.

# Overall Objective:
# To introduce the use of Linear Mixed Models (LMM) for hierarchical data.

# Specific Objectives:
# - Explore the dataset (EDA).
# - Compare linear fixed-effect models with LMM.
# - Evaluate assumptions (normality, homoscedasticity).
# - Interpret results using coefficients, random effects, and R².
# - Visualize predictions and residuals.

# Dataset description:
# growth : Tree growth in diameter (mm/year) - numeric variable
# sp     : Species identifier (e.g., sp01, sp02...) - categorical
# pft    : Plant Functional Type (deciduous / evergreen) - categorical
# height : Tree height (m) - numeric variable

###############################################################
# 1) Libraries (only essential ones)
###############################################################

# List of required packages with explanations
# (Comments are added right after each name)
required_packages <- c(
  "lme4",        # Core package for fitting linear mixed-effects models (LMMs)
  "lmerTest",    # Extends lme4 by providing p-values for fixed effects
  "ggplot2",     # Data visualization
  "car",         # Used here for Levene’s test)
  "performance", # Tools for model evaluation and diagnostics (R², checks, etc.)
  "viridis",     # Colorblind-friendly palettes for plots
  "patchwork"    # Combine multiple ggplot2 plots into one layout
)

# Install any that are missing
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

# Load all libraries
lapply(required_packages, library, character.only = TRUE)


###############################################################
# 2) Import and Explore the Data
###############################################################

# Read the dataset from a CSV file
# sep = ";"   → because the file is semicolon-delimited
# header = TRUE → first row contains column names
dataset <- read.table("U:\\06_Carine\\04_Statistics\\05_Workshop\\script_R\\tree_dataset.csv", 
                      sep = ";", header = TRUE)


# Show the first 6 rows of the dataset to check structure
head(dataset)

# List all unique species codes in the dataset
unique(dataset$sp)

# List the unique plant functional types (deciduous vs evergreen)
unique(dataset$pft)

# Summary statistics of tree growth and height (min, median, mean, max, quartiles)
summary(dataset)

###############################################################
# 2) Exploratory Plots
###############################################################

# Boxplot: Tree growth by plant functional type (PFT)
# Boxes are filled with green (evergreen) and brown (deciduous)
p1 <- ggplot(dataset, aes(x = pft, y = growth, fill = pft)) +
  geom_boxplot(alpha = 0.8, color = "black") +              # black outline, colored fill
  scale_fill_manual(values = c("evergreen" = "darkgreen", 
                               "deciduous" = "saddlebrown")) +  # custom colors
  theme_minimal() +                                         # clean minimal theme
  labs(title = "a) Tree Growth by Plant Functional Type",      # plot title
       x = "Plant Functional Type (PFT)",                   # x-axis label
       y = "Growth (mm/year)")                              # y-axis label

# Boxplot: Distribution of growth by species
# Boxes are filled by plant functional type (pft)
p2<-ggplot(dataset, aes(x = sp, y = growth, fill = pft)) +
  geom_boxplot() +                                          # draw boxplots
  scale_fill_manual(values = c("evergreen" = "darkgreen", 
                               "deciduous" = "saddlebrown")) +  # custom colors
  theme_minimal() +                                         # clean minimal theme
  labs(title = "d) Tree Growth Distribution by Species",       # plot title
       x = "Species", y = "Growth (mm/year)")               # axis labels

# Scatterplot: Tree growth vs. height
# Points are colored by plant functional type (pft)
p3<-ggplot(dataset, aes(x = height, y = growth, color = pft)) +
  geom_point(size = 2, alpha = 0.7) +                       # add points (semi-transparent)
  scale_color_manual(values = c("evergreen" = "darkgreen", 
                               "deciduous" = "saddlebrown")) +  # custom colors
  theme_minimal() +                                         # clean minimal theme
  labs(title = "c) Tree Growth by Height",            # plot title
       x = "Height (m)", y = "Growth (mm/year)")            # axis labels


# Scatterplot: Tree growth vs. height
# Points are colored by species type (sp)
p4 <- ggplot(dataset, aes(x = height, y = growth, color = sp)) +
  geom_point(size = 2, alpha = 0.7) +                          # species points
  theme_minimal() +
  labs(title = "d) Tree Growth by Height (by Species)",         # title
       x = "Height (m)", y = "Growth (mm/year)",                # axis labels
       color = "Species")                                       # legend title


# Combine diagnostic plots in a 2x2 grid
(p1 | p2) / (p3 | p4) #

###############################################################
# 2) Baseline Linear Model (No Random Effects)
###############################################################

# Step 0: Simple linear regression model (fixed effects only)
# This ignores the species grouping structure.
# It assumes all observations are independent, which is not true.
# -> Risk: underestimates standard errors and inflates Type I error.
lm_fixed <- lm(growth ~ height + pft, data = dataset)
summary(lm_fixed)

head(dataset)
dataset$pft<-as.factor(dataset$pft)
model <- lm(growth ~ pft, data = dataset)
summary(model)

library(ggplot2)
ggplot(dataset, aes(x = pft, y = growth)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  ylab("Growth (mm/year)") +
  xlab("Leaf Type")
plot(model)


# Note:
# lm_fixed gives us a baseline, but it does not account for
# species-level variation. We need mixed models to properly handle
# the hierarchical structure (trees nested within species).

###############################################################
# 3) Linear Mixed-Effects Model Selection
###############################################################

# Step 1: Full mixed model (random slope + intercept for species)
# --------------------------------------------------------------
# Fixed effects:
#   - height (continuous predictor of growth, in meters)
#   - pft (plant functional type: deciduous vs evergreen)
#   - height:pft interaction (tests whether the slope of growth vs. height
#     differs between deciduous and evergreen species)
#
# Random effects:
#   - (height | sp): allows each species (sp) to have:
#       * its own intercept (baseline growth)
#       * its own slope with respect to height
#   This structure captures both baseline differences among species
#   and species-specific responses to tree height.
#
# Model formula:
lmm_full <- lmer(growth ~ height + pft + height:pft + (height | sp), data = dataset)
summary(lmm_full)

# --------------------------------------------------------------
# Notes on fixed effects coding:
# --------------------------------------------------------------
# - R encodes categorical variables (factors) with dummy variables.
# - The first level (alphabetical by default) is used as the "reference".
# - Here, pft has two levels: "deciduous" (reference) and "evergreen".
#   → "pftevergreen" means the difference for evergreen species relative
#     to deciduous.
#
# Interpretation of coefficients:
# - (Intercept): expected growth for deciduous species at height = 0.
# - height: slope of growth vs. height for deciduous species.
# - pftevergreen: difference in baseline growth between evergreen
#                 and deciduous species.
# - height:pftevergreen: difference in slope of height effect between
#                        evergreen and deciduous species.


# Step 2: Random intercept model (simpler random structure)
# Fixed effects:
#   - height (continuous predictor of growth)
#   - pft (plant functional type: deciduous vs evergreen)
#   - height:pft interaction (does the effect of height differ by PFT?)
# Random effects:
#   - (1 | sp): only the intercept varies by species,
#     meaning species differ in baseline growth but share the same slope.
# Notes:
#   This is simpler than the random slope model, because it assumes
#   all species respond to height in the same way (parallel lines),
#   but they may start from different baseline levels of growth.
lmm_ri <- lmer(growth ~ height + pft + height:pft + (1 | sp), data = dataset)   # fit LMM with random intercept only
summary(lmm_ri)                                                    # show summary of fixed + random effects

# Compare random slope vs random intercept
anova(lmm_full, lmm_ri, refit = FALSE)                             # likelihood ratio test: slope vs intercept
# If p < 0.05 → keep random slopes
# If p > 0.05 → random intercept is sufficient

# Step 3: Remove interaction (height:pft)
lmm_no_interaction <- lmer(growth ~ height + pft + (height | sp), data = dataset)  # drop interaction term
anova(lmm_full, lmm_no_interaction)                                                # test if interaction is significant
# If p > 0.05 → interaction not needed, keep simpler model

# Step 4: Remove plant functional type effect (pft)
lmm_no_pft <- lmer(growth ~ height + (height | sp), data = dataset)   # drop pft (evergreen/deciduous) effect
anova(lmm_no_pft, lmm_no_interaction)                                 # compare with model including pft
# If p < 0.05 → pft effect matters, keep in model

# Step 5: Remove height effect
lmm_no_height <- lmer(growth ~ pft + (height | sp), data = dataset)   # drop height effect
anova(lmm_no_interaction, lmm_no_height)                              # compare with model including height
# If p < 0.05 → height matters, keep in model

# Final chosen model
final_model <- lmm_no_interaction                                      # assign the selected best model
summary(final_model)                                                   # show results (coefficients, SE, p-values)
fixef(final_model)                                                     # extract fixed effects only

# Explanation:
#   
# Intercept (0.124)
# This is the baseline growth rate (in mm/year) for the reference group (deciduous species) when tree height = 0.
# (Of course, height = 0 is not realistic, but it anchors the model.)
# 
# height (0.0513)
# For every 1 meter increase in tree height, tree growth increases by about 0.051 mm/year (on average, across species).
# This slope applies to both deciduous and evergreen species.
# 
# pftevergreen (0.3623)
# Evergreen species grow, on average, 0.36 mm/year more than deciduous species (the reference group), after accounting for tree height.
# This shows the functional type effect.
# 
# So, the model equation can be written as:
#   
#   Deciduous species:
#   Growth = 0.124 + 0.051 × Height
# 
# Evergreen species:
#   Growth = (0.124 + 0.362) + 0.051 × Height
# Growth = 0.486 + 0.051 × Height

###############################################################
# 4) Model Diagnostics
###############################################################

# Extract residuals (observed - fitted values)
res <- resid(final_model)

# --------------------------------------------------------------
# Shapiro-Wilk normality test
# --------------------------------------------------------------
# Tests whether residuals are normally distributed.
# Hypotheses:
# H0: Residuals follow a normal distribution.
# H1: Residuals deviate from normality.
# If p > 0.05 → assumption of normality is satisfied.
shapiro.test(res)

# --------------------------------------------------------------
# Levene’s test for homogeneity of variance
# --------------------------------------------------------------
# Tests if residual variance is equal across groups (here, PFT).
# Hypotheses:
# H0: Variances are equal across groups (homoscedasticity).
# H1: Variances differ between groups (heteroscedasticity).
# If p > 0.05 → assumption of homogeneity is satisfied.
leveneTest(res ~ as.factor(dataset$pft))

# --------------------------------------------------------------
# R² for mixed model
# --------------------------------------------------------------
# Computes marginal and conditional R² values.
# - Marginal R² = variance explained by fixed effects only.
# - Conditional R² = variance explained by both fixed + random effects.
r2(final_model)

# --------------------------------------------------------------
# Prepare data for plots
# --------------------------------------------------------------
# Build a data frame with fitted values and residuals for visualization.
residual_df <- data.frame(Fitted = fitted(final_model), Residuals = res)

# --------------------------------------------------------------
# Residuals vs Fitted Plot
# --------------------------------------------------------------
# Purpose: check for systematic patterns in residuals.  
# - Ideally: residuals form a random "cloud" around 0 with no trend.  
# - Patterns (curves, funnels) may indicate non-linearity or heteroscedasticity.  
p5 <- ggplot(residual_df, aes(x = Fitted, y = Residuals, color = abs(Residuals))) +  # x = fitted values, y = residuals, color by residual magnitude
  geom_point(size = 2, alpha = 0.8) +                                                # scatterplot of residuals (semi-transparent points)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +                 # dashed reference line at 0 (ideal mean of residuals)
  scale_color_viridis(option = "C", direction = 1) +                                 # viridis scale: blue = small residuals, yellow = large
  theme_minimal() +                                                                  # clean background
  labs(title = "a) Residuals vs Fitted",                                                # plot title
       x = "Fitted Values [Growth (mm/year)]",                                       # x-axis = predicted growth
       y = "Residuals [Growth (mm/year)]",                                           # y-axis = observed - predicted
       color = "|Residuals|")                                                        # legend explains color = |residuals|


# --------------------------------------------------------------
# Histogram of Residuals with Normal Curve
# --------------------------------------------------------------
# Purpose: check approximate normality of residuals.  
# - A bell-shaped histogram suggests residuals are close to normal.  
# - Skewed or multi-modal shapes indicate potential violations.  
p6 <- ggplot(residual_df, aes(x = Residuals)) +  
  geom_histogram(aes(y = after_stat(density), fill = after_stat(x)),                 # histogram of residuals (scaled to density)
                 bins = 20, color = "white", alpha = 0.9) +                          # 20 bins, white borders, semi-transparent fill
  
  stat_function(fun = dnorm,                                                         # overlay theoretical normal density curve
                args = list(mean = mean(res), sd = sd(res)),                         # parameters: mean and sd of residuals
                aes(color = "Theoretical Normal"),                                   # label for legend
                linewidth = 1, linetype = "dashed") +                                # dashed black curve
  
  scale_color_manual(name = "Distribution",                                          # legend for curve
                     values = c("Theoretical Normal" = "black")) +                   # black dashed line
  
  scale_fill_viridis(option = "C", direction = 1) +                                  # viridis color scale for residuals
  
  theme_minimal() +                                                                  # clean background
  
  labs(title = "b) Histogram of Residuals with Normal Curve",                           # plot title
       x = "Residuals [Growth (mm/year)]",                                           # x-axis label
       y = "Density",                                                                # y-axis label
       fill = "Residuals")                                                           # legend for histogram fill


# --------------------------------------------------------------
# Combine diagnostic plots side by side
# --------------------------------------------------------------
(p5 | p6)


###############################################################
# 5) Visualization of Predictions (with equations)
###############################################################

# Predictions from the final model
dataset$pred_full  <- predict(final_model)              # predictions with fixed + random effects (species-specific)
dataset$pred_fixed <- predict(final_model, re.form=NA)  # predictions with fixed effects only (ignores species variation)

# Visualization: observed growth + mixed vs fixed predictions
p4<-ggplot() +
  geom_point(data = dataset,                            # plot observed values
             aes(x = height, y = growth, color = sp),   # x = height, y = growth, color = species
             size = 2, alpha = 0.7) +                   # semi-transparent points
  
  geom_line(data = dataset,                             # plot mixed-effects predictions
            aes(x = height, y = pred_full, group = sp,  # predictions vary by species (grouped)
                color = sp,                             # color lines by species
                linetype = "Mixed (Fixed + Random)"),   # labeled in legend
            size = 0.8) +                               # thinner line for mixed effects
  
  geom_line(data = dataset,                             # plot fixed-effects predictions
            aes(x = height, y = pred_fixed, group = pft,# predictions by plant functional type (pft)
                color = pft,                            # color lines by deciduous/evergreen
                linetype = "Fixed Only"),               # labeled in legend
            size = 1.3) +                               # thicker line for fixed effects
  
  scale_linetype_manual(values = c("Fixed Only" = "solid", 
                                   "Mixed (Fixed + Random)" = "dashed")) + # define line styles
  
  scale_color_manual(                                   # set color palettes
    values = c(
      setNames(turbo(length(unique(dataset$sp))), unique(dataset$sp)), # turbo palette for species
      deciduous = "brown", evergreen = "darkgreen"                     # custom colors for PFT
    )
  ) +
  
  theme_minimal() +                                     # clean theme
  labs(title = "Observed vs Predicted Growth by Species", # plot title
       y = "Growth (mm/year)", x = "Height (m)",         # axis labels
       linetype = "Prediction Type",                     # legend for line types
       color = "Species / PFT") +                        # legend for species and PFT
  
  # Annotate fixed-effect equations directly on the plot
  annotate("text", x = 10,          # position text at right side
           y = 2.4,                # place near bottom
           label = "Deciduous: Growth = 0.124 + 0.051 × Height", # regression equation for deciduous
           color = "brown", hjust = 0, size = 4) +       # brown text for deciduous
  
  annotate("text", x = 10,          # position text at right side
           y = 2.5,                # place near top
           label = "Evergreen: Growth = (0.124 + 0.362) + 0.051 × Height", # regression equation for evergreen
           color = "darkgreen", hjust = 0, size = 4)     # green text for evergreen

print(p4) # plot the graphic 


###############################################################
# END of Workshop Script
# Key Takeaways:
#
# - Linear Mixed-Effects Models (LMMs) handle hierarchical data:
#   Useful when data are grouped (e.g., trees within species, plots, or sites),
#   capturing both fixed effects (population-level) and random effects (group-level).
#
# - Fixed-effect models (LM) ignore grouping structure:
#   They assume all observations are independent, which can underestimate
#   variability and inflate Type I error rates.
#
# - Model selection is based on ANOVA comparisons:
#   Allows testing which fixed and random effects improve the model.
#   Parsimonious models (simpler but informative) are preferred.
#
# - Model diagnostics are essential:
#   Residual plots check linearity and homoscedasticity,
#   Shapiro-Wilk test and histograms assess normality,
#   and Levene’s test checks homogeneity of variances.
#
# - R² interpretation for mixed models:
#   * Marginal R² = variance explained by fixed effects only.
#   * Conditional R² = variance explained by both fixed + random effects.
#   This distinction shows the added contribution of random effects.
#
# - Visualization enhances interpretation:
#   Scatterplots, residual plots, and fitted vs observed predictions
#   clarify how fixed and random effects contribute to explaining variation.
#
# - Practical insight:
#   In this example, both tree height and plant functional type (PFT: deciduous vs evergreen)
#   were informative predictors of tree growth, while random effects captured
#   species-specific deviations around these general trends.
###############################################################
