context("Test cytofQC.")

test_that("basic functionality works", {
    # load small example data
    library(CATALYST)
    data("raw_data")
    
    tech <- readCytof(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
    expect_equal(nrow(tech), 61)
    
    # expect_warning({ lab <- labelQC(tech) },
    #                'Not enough doublets')
    suppressWarnings({ lab <- labelQC(tech) })
    expect_equal(nrow(lab), 61)
    # expect_equal(names(lab), c('Time','label','bead','doublet','debris','dead'))
    # expect_equivalent(table(lab$label)['cell'], 4730)
    expect_equivalent(table(lab$label)['gdpZero'], 270)
    
    
    # load big example data
    f <- system.file("extdata", "raw_cytof.fcs", package = "cytofQC")
    expect_equal(basename(f), "raw_cytof.fcs")
    
    tech <- readCytof(f)
    expect_equal(nrow(tech), 49)
    
    suppressWarnings({ lab <- labelQC(tech, model = 'rf', nTrain = 500) })
    expect_equal(nrow(lab), 49)
    # expect_equal(names(lab), c('Time','label','bead','doublet','debris','dead'))
    expect_equivalent(table(lab$label)['gdpZero'], 1060)
    
})
