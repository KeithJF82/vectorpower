
#### HELPER FUNCTIONS ####################################################################
#------------------------------------------------
#' @title Nice Format
#'
#' @description for single value, return value as string. For vector of values return string of comma-separated values 
#' enclosed in curly brackets
#'
#' @details x is single value or vector of values
#'
#' @param x  Value to convert to string
#'
#' @export
#' 
nice_format <- function(x) {
  if (is.null(x)) {
    return("")
  }
  if (length(x)==1) {
    ret <- as.character(x)
  } else {
    ret <- paste0("{", paste(x, collapse = ", "), "}")
  }
  return(ret)
}

#### BASIC OBJECT TYPES ####################################################################

#------------------------------------------------
#' @title x is not null
#'
#' @description Check that value x is not null
#'
#' @details Check that x is not null
#'
#' @param x        Value to check if not equal to null
#' @param message  Message to output if x is equal to null
#' @param name     Name
#'
#' @export
#' 
assert_non_null <- function(x, message = "%s cannot be null", name = deparse(substitute(x))) {
  if (is.null(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
#' @title x is atomic
#'
#' @description Check that value x is atomic
#'
#' @details Check that x is atomic
#'
#' @param x        Value to check if atomic
#' @param message  Message to output if x not atomic
#' @param name     Name
#'
#' @export
#' 
assert_atomic <- function(x, message = "%s must be atomic (see ?is.atomic)", name = deparse(substitute(x))) {
  if (!is.atomic(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
#' @title x is atomic and single valued (has length 1)
#'
#' @description Check that value x is atomic and single valued (has length 1)
#'
#' @details Check that x is atomic and single valued (has length 1)
#'
#' @param          x Value to check if atomic and single valued (has length 1)
#' @param message  Message to output if x not atomic and single valued
#' @param name     Name
#'
#' @export
#' 
assert_single <- function(x, message = "%s must be a single value", name = deparse(substitute(x))) {
  assert_non_null(x, name = name)
  assert_atomic(x, name = name)
  assert_length(x, 1, name = name)
  return(TRUE)
}

#------------------------------------------------
#' @title Assert String
#'
#' @description x is character string
#'
#' @details x is character string
#'
#' @param x        Value to check if character string
#' @param message  Message to output if x not character string
#' @param name     Name
#'
#' @export
#' 
assert_string <- function(x, message = "%s must be character string", name = deparse(substitute(x))) {
  if (!is.character(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
#' @title Assert String
#'
#' @description x is single character string
#'
#' @details x is single character string
#'
#' @param         x Value to check if single character string
#' @param name    Name
#'
#' @export
#' 
assert_single_string <- function(x, name = deparse(substitute(x))) {
  assert_length(x, n = 1, name = name)
  assert_string(x, name = name)
  return(TRUE)
}

#------------------------------------------------
#' @title Assert logical
#'
#' @description x is logical
#'
#' @details x is logical
#'
#' @param x        Value to check if logical
#' @param message  Message to output if x not logical
#' @param name     Name
#'
#' @export
#' 
assert_logical <- function(x, message = "%s must be logical", name = deparse(substitute(x))) {
  if (!is.logical(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
#' @title Assert Numeric
#'
#' @description x is numeric
#'
#' @details x is numeric
#'
#' @param x        Value to check if numeric
#' @param message  Message to output if x not numeric
#' @param name     Name
#'
#' @export
#' 
assert_numeric <- function(x, message = "%s must be numeric", name = deparse(substitute(x))) {
  if (!is.numeric(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
#' @title Assert Single Numeric
#'
#' @description x is single numeric
#'
#' @details x is single numeric
#'
#' @param x       Value to check if single numeric
#' @param name    Name
#'
#' @export
#' 
assert_single_numeric <- function(x, name = deparse(substitute(x))) {
  assert_length(x, n = 1, name = name)
  assert_numeric(x, name = name)
  return(TRUE)
}

#------------------------------------------------
#' @title Assert Integer
#'
#' @description x is integer
#'
#' @details x is integer
#'
#' @param x        Value to check if integer
#' @param message  Message to output if x not integer
#' @param name     Name
#'
#' @export
#' 
assert_int <- function(x, message = "%s must be integer valued", name = deparse(substitute(x))) {
  assert_numeric(x, name = name)
  if (!isTRUE(all.equal(x, as.integer(x), check.attributes = FALSE))) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
#' @title Assert Single integer
#'
#' @description x is single integer
#'
#' @details x is single integer
#'
#' @param x     Value to check if single integer
#' @param name  Name
#'
#' @export
#' 
assert_single_int <- function(x, name = deparse(substitute(x))) {
  assert_length(x, n = 1, name = name)
  assert_int(x, name = name)
  return(TRUE)
}

#------------------------------------------------
#' @title x is positive (with or without zero allowed)
#'
#' @description Check that value x is positive
#'
#' @details Check that x is positive (zero allowed or disallowed using zero_allowed parameter)
#'
#' @param x             Value to check if positive
#' @param zero_allowed  Indicator (TRUE/FALSE) of whether zero is allowed
#' @param message1      Message to output if x is negative (when zero allowed)
#' @param message2      Message to output if x is not positive (when zero disallowed)
#' @param name          Name
#'
#' @export
#' 
assert_pos <- function(x, zero_allowed = TRUE, message1 = "%s must be greater than or equal to zero", 
                       message2 = "%s must be greater than zero", name = deparse(substitute(x))) {
  assert_numeric(x, name = name)
  if (zero_allowed) {
    if (!all(x>=0)) {
      stop(sprintf(message1, name), call. = FALSE)
    }
  } else {
    if (!all(x>0)) {
      stop(sprintf(message2, name), call. = FALSE)
    }
  }
  return(TRUE)
}

#------------------------------------------------
#' @title x is single positive (with or without zero allowed)
#'
#' @description Check that value x is single positive
#'
#' @details Check that x is single positive (zero allowed or disallowed using zero_allowed parameter)
#'
#' @param x             Value to check if single positive
#' @param zero_allowed  Indicator (TRUE/FALSE) of whether zero is allowed
#' @param name          Name
#'
#' @export
#' 
assert_single_pos <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_length(x, n = 1, name = name)
  assert_pos(x, zero_allowed = zero_allowed, name = name)
  return(TRUE)
}

#------------------------------------------------
#' @title x is positive integer (with or without zero allowed)
#'
#' @description Check that value x is positive integer
#'
#' @details Check that x is positive integer (zero allowed or disallowed using zero_allowed parameter)
#'
#' @param x            Value to check if positive integer
#' @param zero_allowed Indicator (TRUE/FALSE) of whether zero is allowed
#' @param name         Name
#'
#' @export
#' 
assert_pos_int <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_int(x, name = name)
  assert_pos(x, zero_allowed = zero_allowed, name = name)
  return(TRUE)
}

#------------------------------------------------
#' @title x is single positive integer (with or without zero allowed)
#'
#' @description Check that value x is single positive integer
#'
#' @details Check that x is single positive integer (zero allowed or disallowed using zero_allowed parameter)
#'
#' @param x             Value to check if single positive integer
#' @param zero_allowed  Indicator (TRUE/FALSE) of whether zero is allowed
#' @param name          Name
#'
#' @export
#' 
assert_single_pos_int <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_length(x, n = 1, name = name)
  assert_pos_int(x, zero_allowed = zero_allowed, name = name)
  return(TRUE)
}

#------------------------------------------------
#' @title Assert Single bounded
#'
#' @description x is single value bounded between limits
#'
#' @details x is single value bounded between limits
#'
#' @param x                 value to check if single value bounded between limits
#' @param left              left bound
#' @param right             right bound
#' @param inclusive_left    True/False is left bound inclusive?
#' @param inclusive_right   True/False is left bound inclusive?
#' @param name              Name
#'
#' @export
#' 
assert_single_bounded <- function(x, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE, 
                                  name = deparse(substitute(x))) {
  assert_length(x, n = 1, name = name)
  assert_bounded(x, left = left, right = right, inclusive_left = inclusive_left, inclusive_right = inclusive_right, name = name)
  return(TRUE)
}

#------------------------------------------------
#' @title Assert Vector
#'
#' @description x is a vector (and is not a list or another recursive type)
#'
#' @details x is a vector (and is not a list or another recursive type)
#'
#' @param x       Value to check if vector and not recursive
#' @param message Message to output if x not a non-recursive vector
#' @param name    Name
#'
#' @export
assert_vector <- function(x, message = "%s must be a non-recursive vector", name = deparse(substitute(x))) {
  if (!is.vector(x) || is.recursive(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
#' @title Assert list
#'
#' @description x is a list
#'
#' @details x is a list
#'
#' @param x       Value to check if list
#' @param message Message to output if x not a list
#' @param name    Name
#'
#' @export
assert_list <- function(x, message = "%s must be a list", name = deparse(substitute(x))) {
  if (!is.list(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}


#### VALUE COMPARISONS ####################################################################

#------------------------------------------------
# x is between bounds (inclusive or exclusive)
#' @title Assert that x is between bounds
#'
#' @description x is between bounds
#'
#' @details x is between bounds (inclusive or exclusive)
#'
#' @param x               Value to verify as being within bounds
#' @param left            Lower bound
#' @param right           Upper bound
#' @param inclusive_left  True/false indicator if lower bound is inclusive
#' @param inclusive_right True/false indicator if upper bound is inclusive
#' @param name            Name
#'
#' @export
assert_bounded <- function(x, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE, 
                           name = deparse(substitute(x))) {
  assert_numeric(x, name = name)
  if (inclusive_left) {
    if (!all(x>=left)) {
      stop(sprintf("%s must be greater than or equal to %s", name, left), call. = FALSE)
    }
  } else {
    if (!all(x>left)) {
      stop(sprintf("%s must be greater than %s", name, left), call. = FALSE)
    }
  }
  if (inclusive_right) {
    if (!all(x<=right)) {
      stop(sprintf("%s must be less than or equal to %s", name, right), call. = FALSE)
    }
  } else {
    if (!all(x<right)) {
      stop(sprintf("%s must be less than %s", name, right), call. = FALSE)
    }
  }
  return(TRUE)
}

#------------------------------------------------
#' @title Assert that x is a subset of y
#'
#' @description x is a subset of y
#'
#' @details x is a subset of y (all x are in y)
#'
#' @param x         Value(s) to check if subset of y
#' @param y         Set of values x must be a subset of
#' @param message   Message to output if not all x are in y
#' @param name_x    Name of x
#' @param name_y    Name of y
#'
#' @export
assert_in <- function(x, y, message = "all %s must be in %s",
                      name_x = deparse(substitute(x)), name_y = nice_format(y)) {
  assert_non_null(x, name = name_x)
  assert_non_null(y, name = name_y)
  if (!all(x %in% y)) {
    stop(sprintf(message, name_x, name_y), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# none of x are in y
#' @title Assert that none of x are in y
#'
#' @description None of x are in y
#'
#' @details None of x are in y
#'
#' @param x         Value(s) to check are not in y
#' @param y         Set of values x must not be in
#' @param message   Message to output if some of x are in y
#' @param name_x    Name of x
#' @param name_y    Name of y
#'
#' @export
assert_not_in <- function(x, y, message = "none of %s can be in %s",
                          name_x = deparse(substitute(x)), name_y = nice_format(y)) {
  assert_non_null(x, name = name_x)
  assert_non_null(y, name = name_y)
  if (any(x %in% y)) {
    stop(sprintf(message, name_x, name_y), call. = FALSE)
  }
  return(TRUE)
}


#### DIMENSIONS ####################################################################

#------------------------------------------------
# length(x) equals n
#' @title Assert that x is of a specified length n
#'
#' @description x is of a specified length n
#'
#' @details x is of a specified length n
#'
#' @param x         Vector to check length of
#' @param n         Required length of x
#' @param message   Message to output if x is not of length n
#' @param name      Name
#'
#' @export
assert_length <- function(x, n, message = "%s must be of length %s", name = deparse(substitute(x))) {
  assert_pos_int(n)
  if (length(x) != n[1]) {
    stop(sprintf(message, name, n), call. = FALSE)
  }
  return(TRUE)
}

#### MISC ####################################################################

#------------------------------------------------
#' @title Check that file exists at chosen path
#'
#' @description Check that file exists at chosen path
#'
#' @details Check that file exists at chosen path
#'
#' @param x         File location
#' @param message   Message to output if file does not exist
#' @param name      Name
#'
#' @export
#' 
assert_file_exists <- function(x, message = "file not found at path %s", name = deparse(substitute(x))) {
  if (!file.exists(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
#' @title Find position
#'
#' @description Find position of value in set
#'
#' @details Find best position of value in set by finding closest value; modification of FindInterval() function
#'
#' @param value   Value to locate
#' @param set     Set of values in which to locate value
#'
#' @export
#' 
findPosition <- function(value=0,set=c(0,1)) {
  j=max(1,findInterval(value,set))
  if(j<length(set)){
    value_a=set[j]
    value_b=set[j+1]
    if(value-value_a>value_b-value){j=j+1}
  }
  return(j)
}
