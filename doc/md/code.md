# Code Specification

References:
- [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html)
- [Doxygen documentation](https://www.doxygen.nl/manual/index.html)

## Header Files

1. Every `.cpp` file should have an associated `.h` file, except unit tests and `main.cpp` containing just a `main()` function.
2. All header files should have `#define` guards to prevent multiple inclusion. The format of the macro should be `<PROJECT>_<PATH>_<FILE>_H_`.
3. All of a project's header files should be listed as descendants of the project's source directory.
4. Order your includes as follows:
   - Associated header of the file.
   - C system headers.
   - C++ standard library headers.
   - Other libraries' .h files.
   - Your project's .h files.
5. Namespaces subdivide the global scope into distinct, named scopes, and so are useful for preventing name collisions in the global scope. Give each file a unique namespace (its path).
6. Use unnamed namespace in `.cpp` for functions which are not contained in associated header file instead of `static` (or even do nothing).

## Variable Names

1. The names of namespace, class, struct and enum should start with a capital letter and have a capital letter for each new word, with no underscores, such as `MyExcitingClass`.
2. The names of variables, functions and data members are all lowercase, with underscores between words, such as `a_struct_data_member`.
3. Data members of classes have additional underscores according to its access permissions. For instance: `__private_member`, `_protected_member`, `public_member`. Don't use public member as possible as you can. Don't use `this->` when access members except members of parent class.
4. The names of macro and constant are with all capitals and underscores, such as `PI_CONST`.

## Comments and Formats

1. Comments are in doxygen format, such as:
   - `/// something` for single line comment.
   - Every file and class (in header file) should have its corresponding description, use `@name`.
2. 