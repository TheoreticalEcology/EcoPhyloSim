// Shared script of the EcoPhyloSim vignettes, added by a js chunk: shows the abstract
// before the table of contents (rmarkdown::html_vignette places the contents first)
document.addEventListener("DOMContentLoaded", function () {
  var toc = document.getElementById("TOC"), abstract = document.querySelector(".abstract");
  if (toc && abstract) abstract.after(toc);
});
