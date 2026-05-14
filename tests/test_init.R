
library(hwas)

source(system.file("utils/unittest.R", package = "hwas"))

unittests$TEST("InitProject",
               "test_is_valid_dirname", 
               function() {
    Expect$true( hwas::is_valid_dirname("grms"))
    Expect$true( hwas::is_valid_dirname("trait/values"))
    Expect$true( hwas::is_valid_dirname("trait_values"))
    Expect$true( hwas::is_valid_dirname("trait-values"))
    Expect$false(hwas::is_valid_dirname("trait.values"))
    Expect$false(hwas::is_valid_dirname("trait?values"))
    Expect$false(hwas::is_valid_dirname("trait\values"))
})


if (!interactive())
    unittests$main()
