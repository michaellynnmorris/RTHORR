context("RTHORR testing")

testthat::test_that("RTHORR testing",

                    {

                      # previous results
                      previous_randall_output <- readr::read_rds(system.file("extdata", "randall_output_test.rds", package = "RTHORR"))
                      previous_randmf_output <- readr::read_rds(system.file("extdata", "randmf_output_test.rds", package = "RTHORR"))

                      previous_randall_from_df_output <- readr::read_rds(system.file("extdata", "previous_randall_from_df_output.rds", package = "RTHORR"))
                      previous_randmf_from_df_output <- readr::read_rds(system.file("extdata", "previous_randmf_from_df_output.rds", package = "RTHORR"))


                      #read correlation matrixes for testing
                      input_matrixes <- system.file("extdata", "input.txt", package = "RTHORR")

                      #read random dataframes for testing
                      df_list <- readr::read_rds(system.file("extdata", "df_list_for_testing.rds", package = "RTHORR"))



                      #run randall on input.txt
                      #outputs a single data frame with RTHOR results
                      randall_output <- RTHORR::randall(n=6,
                                                        nmat=3,
                                                        ord = "circular6",
                                                        input=input_matrixes,
                                                        description = c("sample_one", "sample_two", "sample_three"))


                      #run randmf on input.txt
                      #outputs a list with two dataframes (RTHOR results and comparisons)
                      randmf_output <- RTHORR::randmf(n=6,
                                                      nmat=3,
                                                      ord = "circular6",
                                                      input=input_matrixes)



                      #run randall_from_df on raw data
                      test <- randall_from_df(df_list = df_list, description = c("whole sample",
                                                                                 "t1",
                                                                                 "t2",
                                                                                 "t3",
                                                                                 "t4"))


                      # run randmf_from_df on raw data
                      test3 <- randmf_from_df(df_list = df_list)




                      #compare previous results against new results
                      testthat::expect_identical(previous_randmf_output, randmf_output)
                      testthat::expect_identical(previous_randall_output, randall_output)
                      testthat::expect_identical(test, previous_randall_from_df_output)
                      testthat::expect_identical(test3, previous_randmf_from_df_output)


                    }
)
