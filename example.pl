# Call `polymake --script example.pl` after loading the extension, to run this example.

use application "common";
declare $matrix = new Matrix<Integer> (<<".");
 0 1 1 1 1 1
 1 1 0 0 1 1
 1 1 1 0 0 1
 1 0 1 1 0 1
 1 0 0 1 1 1
.
declare $rows = new Vector<Integer> ();
declare $cols = new Vector<Integer> ();
declare $result = is_totally_unimodular ($matrix, $rows, $cols);
if (is_totally_unimodular($matrix, $rows, $cols))
{
  print "The matrix is totally unimodular.\n";
}
else
{
  print "The matrix is NOT totally unimodular.\n";
  print "Rows of the violating submatrix:    $rows\n";
  print "Columns of the violating submatrix: $cols\n";
}
