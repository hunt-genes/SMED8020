# Module 0: Introduction to the Ubuntu Terminal and R

## Learning objectives

By the end of this module, you will be able to:

- Navigate the Ubuntu Terminal (command line interface)
- Execute command-line programs including **bash functions**, **PLINK**, and **R**
- Understand how R fits into a terminal-based data analysis workflow
- Run a basic R script from the Terminal using input and output files

---

## Part 1: Introduction to the Ubuntu Terminal

### What is the Terminal?

HUNT Cloud uses the Linux operating system **Ubuntu**. One of the primary ways to interact with Linux is through the **Terminal**, which provides a **command line interface (CLI)**.

In a CLI:
- You type commands
- The shell interprets them
- The operating system executes them and returns output

The most common Linux shell is **bash** (Bourne Again SHell).

---

### Opening the Terminal

1. Go to the **Launcher**
2. Select **Terminal**

You should now see a prompt similar to:

```bash
username@machine-name:~$
```

---

### Basic Terminal commands

Try running the following commands:

```bash
ls          # list files
pwd         # print working directory
mkdir test  # create a directory
cd test     # move into the directory
cd ..       # move back up one level
```

---

## Part 2: Executing Commands and Programs from the Terminal

One of the most powerful features of the Terminal is the ability to **execute programs**.

If a command is available on your system, you can run it by typing its name.

---

### Example 1: A simple Bash function

Define a bash function:

```bash
hello () {
  echo "Hello from bash!"
}
```

Run the function:

```bash
hello
```

---

### Example 2: Running PLINK from the Terminal

Check that PLINK is available:

```bash
plink --help
```

Example PLINK command:

```bash
plink   --bfile mydata   --freq   --out mydata_freq
```

---

### Example 3: Running R from the Terminal

Start an interactive R session:

```bash
R
```

Exit R:

```r
q()
```

---

## Part 3: Running an R Script from the Terminal

R scripts can take **command-line input** and write output files.

### Example: `hello.R`

Create a file called `hello.R` with the following contents:

```r
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript hello.R <name> <output_file>")
}

name <- args[1]
outfile <- args[2]

message <- paste("Hello", name, "from R!")

cat(message, file = outfile)
cat(message, "
")
```

Make the script executable (optional):

```bash
chmod +x hello.R
```

Run the script from the Terminal:

```bash
Rscript hello.R Alice hello_output.txt
```

View the output:

```bash
cat hello_output.txt
```

---

## Summary

In this module you learned:

- How to navigate the Ubuntu Terminal
- How bash functions, PLINK, and R are all executed from the command line
- How to run R both interactively and via scripts
- How to pass input to an R script and produce output files

This foundation will support reproducible and scalable data analysis workflows.
