library('caret')
library('mlr')

# Load raw data and remove debugging metadata
raw_data = read.csv("../data/subg_ml_with_metadata.csv")
data = raw_data[, !names(raw_data) %in% c("Project", "Accession")]

# Reformat target variable for the mlr library (workaround for bug in mlr)
data$Disease = as.character(data$Disease)

# Define objective; classify datapoints into classes defined by variable Disease
task = makeClassifTask(id="diseased", data=data, target="Disease", positive = "perid")

# Define learner algorithms with optimised hyperparameters
lrns = list(
    makeLearner(id="randomForest", "classif.randomForest", mtry=30, nodesize=11, maxnodes=10000, predict.type='prob'),
    makeLearner(id="boosting", "classif.boosting", coeflearn='Zhu', mfinal=53, predict.type='prob'),
    makeLearner(id="cforest", "classif.cforest", ntree=116, mtry=30, mincriterion=0.889, minbucket=1, fraction=0.298, predict.type='prob'),
    makeLearner(id="ranger", "classif.ranger", mtry=2014,  min.node.size=5, predict.type='prob')
    )

# Run cross-validation using stratified sets of samples
rdesc = makeResampleDesc("CV", iters = 5L, stratify=TRUE)
# Define performance metrics
meas = list(acc, ppv, npv)
# Run benchmark
bmr = benchmark(lrns, task, rdesc, measures=meas)

# Retrieve benchmark performances
bench = getBMRPerformances(bmr, as.df=TRUE)

# Box plots for benchmarked algorithms
bench$learner = bench$learner.id

# Boxplot of algorithm performance
ggplot(bench, aes(x=learner, y=acc, colour=learner)) + geom_boxplot() + ylim(0,1)

# Generate ROC curves
roc_r = generateThreshVsPerfData(bmr, list(fpr, tpr, auc), aggregate = FALSE)
roc_r$data$iter <- as.factor(roc_r$data$iter)
plotROCCurves(roc_r) + aes(color=learner) + facet_grid(~ learner)
