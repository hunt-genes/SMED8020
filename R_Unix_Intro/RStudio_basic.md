# Working in RStudio: A Guided Introduction to R

This section is written for students to follow **step by step inside RStudio**. You will type commands into the **R Console** or a **script file** and observe the output.

---

## 1. Getting started in RStudio

Open RStudio. You should see:
- The **Console** (where R runs commands)
- One or more editor panes (where scripts can be written)

You can type commands directly into the Console, but for reproducibility it is best to **create an R script**:

**File → New File → R Script**

---

## 2. Defining variables and data types

R has five basic (atomic) classes of objects:

- `character`
- `numeric` (real numbers)
- `integer`
- `complex`
- `logical` (`TRUE` / `FALSE`)

### Assigning variables

Use the arrow `(<-)` to assign values:

```r
x <- 5
```

Printing values:

```r
x            # automatic printing
print(x)     # explicit printing
```

Check the class of an object:

```r
class(x)
```

### Coercing between types

```r
y <- as.integer(x)
class(y)

z <- 1L     # L forces integer type
class(z)
```

---

## 3. Vectors

The most basic object in R is a **vector**. Vectors can contain only one data type.

### Creating vectors

```r
# initialize a numeric vector of length 10
a <- vector("numeric", 10)
print(a)
```

```r
# create a vector using c()
b <- c(1, 2, 3)
print(b)
```

Check their classes:

```r
class(a)
class(b)
```

---

## 4. Data frames

A **data frame** stores tabular data.

```r
x <- 1:10
y <- rnorm(10)

df <- data.frame(x, y)
```

Explore the data frame:

```r
dim(df)      # dimensions
names(df)    # column names
```

---

## 5. Factors (categorical data)

Factors represent categorical variables and are common in statistical models.

```r
x <- factor(c("male", "female", "male", "female", "female", "female"))
```

Inspect the factor:

```r
factor(x)
table(x)
```

---

## 6. Missing values

Missing data in R are represented by `NA` or `NaN`.

```r
y <- factor(c(x, NA))
class(y)
```

Identify missing values:

```r
is.na(y)
```

Count values, including missing:

```r
table(y, useNA = "always")
```

Handling missing values in calculations:

```r
z <- c(1, 20, 13, NA, 45)
mean(z)
mean(z, na.rm = TRUE)
```

---

## 7. Libraries and packages

Packages extend R with additional functions, datasets, and tools.

Because of the HUNT Cloud setup, packages are installed using **conda** in the Terminal, not from inside R.

Once installed, packages are loaded in R using `library()`.

Check your R session information:

```r
sessionInfo()
```

Load the `ggplot2` package:

```r
library(ggplot2)
```

---

## 8. Using functions from a package

Create a simple scatter plot using `ggplot2`:

```r
ggplot(df, aes(x, y)) +
  geom_point()
```

This command:
- Uses `df` as the data source
- Maps `x` and `y` to the plot axes
- Draws points using `geom_point()`

---

## Key takeaway

In this section you learned how to:

- Work interactively in RStudio
- Create and inspect variables
- Use vectors, data frames, and factors
- Handle missing data
- Load packages and use functions from them

These fundamentals are essential for all later data analysis work in R.

[Jump back to the main tutorial](https://github.com/hunt-genes/SMED8020/blob/main/R_Unix_Intro/README.md#part-4-running-an-r-script-from-the-terminal)
