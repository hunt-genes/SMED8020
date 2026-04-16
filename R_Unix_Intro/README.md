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

## Part 2: Working with Files and Text in the Terminal

This section demonstrates how to inspect, search, copy, rename, and delete files using a sample text file called `test.txt`.

### Contents of `test.txt`

The example file contains the following text:

```text
Hello world.
Ready to learn R?
Now you are using Terminal!
```

---

## Viewing Files

### Using `less`

```bash
less test.txt
```

- Opens the file in an interactive viewer
- Allows scrolling without printing everything at once
- Press `q` to quit

Use `less` for safely inspecting large files.

---

### Using `cat`

```bash
cat test.txt
```

- Prints the entire contents of the file to the screen
- Best for small files

---

## Searching Files with `grep`

```bash
grep "R" test.txt
```

Output:

```text
Ready to learn R?
```

- Searches for lines containing the letter `R`
- `grep` is case-sensitive by default

Useful options:

```bash
grep -i "r" test.txt   # ignore case
grep -n "R" test.txt   # show line numbers
```

---

## Printing Text with `echo`

```bash
echo "5"
```

- Prints text to the Terminal
- Commonly used in scripts and for testing commands

You can also write text to a file:

```bash
echo "5" > number.txt
```

---

## Copying, Renaming, and Deleting Files

### Copying a File with `cp`

```bash
cp test.txt test2.txt
```

- Creates a copy of `test.txt` called `test2.txt`
- The original file is unchanged

---

### Renaming a File with `mv`

```bash
mv test2.txt test3.txt
```

- Renames `test2.txt` to `test3.txt`
- Also used to move files between directories

---

### Deleting a File with `rm`

```bash
rm test3.txt
```

- Permanently deletes the file
- There is no undo

Safer option:

```bash
rm -i test3.txt
```

---

## Piping Commands Together

### What is a Pipe?

The pipe operator `|` sends the **output of one command** directly into another command.

This allows you to build powerful command chains without creating intermediate files.

---

### Example: `cat` + `grep`

```bash
cat test.txt | grep "R"
```

Explanation:

1. `cat test.txt` outputs the file contents
2. `|` passes that output to `grep`
3. `grep "R"` filters lines containing `R`

Output:

```text
Ready to learn R?
```

This produces the same result as:

```bash
grep "R" test.txt
```

But piping becomes especially useful when combining many commands.

---

### Example: Counting Matching Lines

```bash
grep "R" test.txt | wc -l
```

- `grep "R" test.txt` finds matching lines
- `wc -l` counts how many lines were matched

Output:

```text
1
```

---

## Key Takeaway

Linux commands are small tools designed to work together:

- Files can be inspected with `less` and `cat`
- Content can be searched with `grep`
- Files can be copied, renamed, and deleted with `cp`, `mv`, and `rm`
- Pipes (`|`) allow commands to be chained into powerful workflows

These concepts are essential for working with data, logs, and outputs from tools like PLINK and R.   

## Part 3: Executing Commands and Programs from the Terminal

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

## Part 4: Running an R Script from the Terminal

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
