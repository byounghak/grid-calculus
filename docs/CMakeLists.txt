# docs/CMakeLists.txt
#
# Documentation build pipeline. Produces three artifacts:
#   - Doxygen (HTML + LaTeX) under build/docs/doxygen/
#   - User Guide PDF       under build/docs/user-guide/user-guide.pdf
#   - Developer Note PDF   under build/docs/developer-note/developer-note.pdf
#
# The Developer Note's LaTeX source `\input{}`s the Doxygen-generated
# refman.tex from build/docs/doxygen/latex/, so the developer-note target
# depends on the doxygen target. Each book's source tree is copied into
# build/docs/<book>/ at configure time so that latexmk operates on a clean
# tree alongside the Doxygen output it inputs.

find_package(Doxygen REQUIRED)
find_program(LATEXMK_EXECUTABLE latexmk REQUIRED)

# Configure the Doxyfile from the template, substituting CMake variables.
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in"
    "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
    @ONLY
)

# -- gridcalc-docs-doxygen -----------------------------------------------------

add_custom_target(gridcalc-docs-doxygen
    COMMAND ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Running Doxygen over public headers"
    VERBATIM
)

# -- gridcalc-docs-userguide ---------------------------------------------------

# Copy User Guide source tree into build/docs/user-guide/ so latexmk can
# write its build artifacts alongside the (read-only) source.
add_custom_target(gridcalc-docs-userguide-prep
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/user-guide"
        "${CMAKE_BINARY_DIR}/docs/user-guide"
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/common"
        "${CMAKE_BINARY_DIR}/docs/common"
    COMMENT "Staging User Guide source tree under build/docs/user-guide/"
    VERBATIM
)

# latexmk exits non-zero for benign warnings (undefined refs on first pass,
# Doxygen-generated multiply-defined labels). Wrap in a shell command that
# treats "PDF was produced and is non-empty" as success.
add_custom_target(gridcalc-docs-userguide
    COMMAND sh -c "${LATEXMK_EXECUTABLE} -xelatex -interaction=nonstopmode -f -file-line-error main.tex || true && test -s main.pdf"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/docs/user-guide"
    DEPENDS gridcalc-docs-userguide-prep
    COMMENT "Building User Guide PDF (latexmk)"
    VERBATIM
)

# -- gridcalc-docs-developernote ----------------------------------------------

add_custom_target(gridcalc-docs-developernote-prep
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/developer-note"
        "${CMAKE_BINARY_DIR}/docs/developer-note"
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/common"
        "${CMAKE_BINARY_DIR}/docs/common"
    COMMENT "Staging Developer Note source tree under build/docs/developer-note/"
    VERBATIM
)

# Extract the API-reference content (the per-class \input{} list) from
# Doxygen's standalone refman.tex into a chapter-only fragment that the
# Developer Note can \input{} from inside its own \appendix section.
add_custom_target(gridcalc-docs-extract-api-ref
    COMMAND ${CMAKE_COMMAND} -E env
        sh "${CMAKE_CURRENT_SOURCE_DIR}/extract-api-reference.sh"
        "${CMAKE_BINARY_DIR}/docs/doxygen/latex/refman.tex"
        "${CMAKE_BINARY_DIR}/docs/developer-note/api-reference.tex"
    DEPENDS gridcalc-docs-doxygen gridcalc-docs-developernote-prep
    COMMENT "Extracting API-reference content from refman.tex"
    VERBATIM
)

# TEXINPUTS prepends Doxygen's latex/ output to the LaTeX search path so
# that \input{classgridcalc_...} from api-reference.tex resolves to the
# Doxygen-generated per-class .tex files. The trailing `:` keeps the
# default paths. The `|| true && test -s main.pdf` idiom handles latexmk's
# non-zero exit on benign warnings (multiply-defined labels in the
# Doxygen output, undefined refs on the first pass before convergence).
add_custom_target(gridcalc-docs-developernote
    COMMAND sh -c "TEXINPUTS=\"${CMAKE_BINARY_DIR}/docs/doxygen/latex:\" ${LATEXMK_EXECUTABLE} -xelatex -interaction=nonstopmode -f -file-line-error main.tex || true && test -s main.pdf"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/docs/developer-note"
    DEPENDS gridcalc-docs-developernote-prep
            gridcalc-docs-doxygen
            gridcalc-docs-extract-api-ref
    COMMENT "Building Developer Note PDF (latexmk; inputs Doxygen refman content)"
    VERBATIM
)

# -- gridcalc-docs-all (umbrella) ---------------------------------------------

add_custom_target(gridcalc-docs-all
    DEPENDS gridcalc-docs-doxygen
            gridcalc-docs-userguide
            gridcalc-docs-developernote
    COMMENT "Building all gridcalc documentation artifacts"
)
