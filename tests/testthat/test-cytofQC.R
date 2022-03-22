context("Test cytofQC.")

test_that("basic functionality works", {
    # load example data
    library(CATALYST)
    data("raw_data")
    
    tech <- dataPrep(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
    expect_equal(nrow(tech), 5000)
    
    lab <- labelQC(tech)
    expect_equal(nrow(lab), 5000)
    expect_equal(names(lab), c('Time','label','bead','doublet','debris','dead'))

})
