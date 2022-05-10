<html>



<head>
<title>Diagrams of <?php echo $_GET["knot"] ?>  </title>
</head>

<center>  <h3> Diagrams of <?php echo $_GET["knot"] ?> </h3> </center>

<hr>

<table border="0" cellspacing = "30">

<tr height = "10">
<td><img width="400" height="400" align="center" src="diagrams/<?php echo $_GET["knot"] ?>.png"></td>
<td><img width="400" height="400" align="center" src="diagrams/<?php echo $_GET["knot"] ?>mirror.png"></td>
</tr>

<tr>
<td align="center"><?php echo $_GET["knot"] ?> </td>
<td align="center"><?php echo $_GET["knot"] ?>  mirror</td>
</tr>


</table>
