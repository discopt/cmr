# Call `polymake --script example.pl` after loading the extension, to run this example.

use application "matroid";
declare $matrix = new Matrix<Integer> (<<".");
 1 1 1 1
 1 1 0 0
 1 0 1 0
 1 0 0 1
.
declare $rows = new Vector<Integer> ();
declare $cols = new Vector<Integer> ();
if (is_totally_unimodular($matrix, $rows, $cols))
{
  print "The matrix is totally unimodular.\n";
}
else
{
  print "The matrix is NOT totally unimodular.\n";
  print "  Rows of the violating submatrix:    $rows\n";
  print "  Columns of the violating submatrix: $cols\n";
}
if (is_unimodular($matrix))
{
  print "The matrix is unimodular.\n";
}
else
{
  print "The matrix is NOT unimodular.\n";
}
if (is_strongly_unimodular($matrix))
{
  print "The matrix is strongly unimodular.\n";
}
else
{
  print "The matrix is NOT strongly unimodular.\n";
}
declare $k = is_k_modular($matrix);
if ($k)
{
  print "The matrix is $k-modular.\n";
}
else
{
  print "The matrix is NOT k-modular.\n";
}
if ($k = is_strongly_k_modular($matrix))
{
  print "The matrix is strongly $k-modular.\n";
}
else
{
  print "The matrix is NOT strongly k-modular.\n";
}
