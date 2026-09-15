library(epiworldR)


#' Factory function for wastewater detection
#' and intervention
#' @param m Model to which the intervention will be added.
#' @param mean_shedding_case Mean shedding per infected case.
#' @param detect_prob Detection probability per shed unit.
#' @param intervention_threshold Threshold at which the detected
#' shedding triggers a reduction in contact
#' @param crate_param_name Name of the contact rate param of the
#' model.
#' @param crate_reduction Contact rate reduction rate.
#' @param sampling_interval Interval at which the wastewater is sampled.
wastewater_detect_factory <- function(
  m,
  mean_shedding_case     = 10,
  detect_prob            = .1,
  intervention_threshold = 20,
  crate_param_name       = "Contact rate",
  crate_reduction        = .2,
  sampling_interval      = 7
) {

  crate_param_name_base  <- "Contact rate (base) WW"
  crate_param_name_redux <- "Contact rate (reduced) WW"

  # Adding the corresponding parameters
  add_param(
    m,
    pname = crate_param_name_base,
    pval = get_param(
      m,
      crate_param_name
    )
  )

  add_param(
    m,
    pname = crate_param_name_redux,
    pval = get_param(
      m,
      crate_param_name
    ) * (1 - crate_reduction)
  )
  
  add_param(
    m,
    pname = "Detection is active",
    pval = -1
  )

  fun <- function(m) {

    # Setting the active status to -1 and the
    # original contact rate if it is the first day
    day <- today(m)
    if (day == 1) {
      set_param(
        x = m,
        pname = crate_param_name,
        get_param(m, crate_param_name_base)
      )
      set_param(
        x = m,
        pname = "Detection is active",
        pval = -1
      )

      # Creating the tmp output
      assign("wastewater_counter", value = c(), envir = .GlobalEnv)
    }

    # Should we sample today?
    if (day %% sampling_interval != 0)
      return(NULL)

    # Get status
    is_active <- get_param(m, "Detection is active") > 0

    # Extracting the number of cases historically
    today_infected <- get_today_total(m)["Infected"]

    # The function always returns void
    if (today_infected < 1)
      return(NULL)

    # Using the number of infected as a parameter of
    detected <- rbinom(
      n = 1,
      size = today_infected * mean_shedding_case,
      prob = detect_prob
      ) / mean_shedding_case

    # If the number of cases is above, then
    # we induce an intervention that reduces the contact
    # rate by 20%
    if (!is_active && (detected > intervention_threshold)) {
      set_param(
        x = m,
        pname = crate_param_name,
        get_param(m, crate_param_name_redux)
        )
      set_param(
        x = m,
        pname = "Detection is active",
        pval = 1
      )

      assign(
        "wastewater_counter",
        value = c(
          get("wastewater_counter", envir = .GlobalEnv),
          today(m)
          ),
        envir = .GlobalEnv
      )

    } else if (is_active && (detected <= intervention_threshold)) {
      set_param(
        x = m,
        pname = crate_param_name,
        get_param(m, crate_param_name_base)
        )
      set_param(
        x = m,
        pname = "Detection is active",
        pval = -1
      )

      assign(
        "wastewater_counter",
        value = c(
          get("wastewater_counter", envir = .GlobalEnv),
          -today(m)
          ),
        envir = .GlobalEnv
      )
    }
  

  } 
  
  event <- globalevent_fun(
    fun = fun,
    name = "Wastewater detection",
    day = -99
  )

  add_globalevent(
    m,
    event
  )

  invisible(m)

}
