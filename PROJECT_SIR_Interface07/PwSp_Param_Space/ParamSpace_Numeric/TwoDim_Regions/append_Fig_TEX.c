#include "../../../SIR.h"
#include "../../../SIR_Analytic_General.h"

void append_Fig_TEX(FILE *F_TEX, char *name, float imm)
{

  fprintf(F_TEX, "\n\\begin{figure}");
  fprintf(F_TEX, "\n\\centering");
  fprintf(F_TEX, 
	  "\n{\\par\\centering \\resizebox*{0.7\\columnwidth}{!}{\\rotatebox{0}");
  fprintf(F_TEX, "\n{\\includegraphics{../%s}}}\\par}", name);
  fprintf(F_TEX, "\n\\caption{File: %s. Regions in parameter Space, $\\eta = %f$}", name, imm);
  fprintf(F_TEX, "\n\\end{figure}");
  fprintf(F_TEX, "\n\n\\pagebreak");
}

