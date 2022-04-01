context("Test cytofQC.")

test_that("basic functionality works", {
    # load small example data
    library(CATALYST)
    data("raw_data")
    
    tech <- dataPrep(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
    expect_equal(nrow(tech), 5000)
    
    expect_warning({ lab <- labelQC(tech) },
                   'Not enough doublets')
    expect_equal(nrow(lab), 5000)
    expect_equal(names(lab), c('Time','label','bead','doublet','debris','dead'))
    expect_equivalent(table(lab$label)['cell'], 4730)
    expect_equivalent(table(lab$label)['GDPzero'], 270)
    
    
    # load big example data
    f <- system.file("extdata", "raw_cytof.fcs", package = "cytofQC")
    expect_equal(basename(f), "raw_cytof.fcs")
    
    tech <- dataPrep(f)
    expect_equal(nrow(tech), 30000)
    
    expect_warning({ lab <- labelQC(tech, model = 'rf', nTrain = 500) },
                   'Not enough dead')
    expect_equal(nrow(lab), 30000)
    expect_equal(names(lab), c('Time','label','bead','doublet','debris','dead'))
    expect_equivalent(table(lab$label)['GDPzero'], 1060)
    
})
