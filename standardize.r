# standardize
standardize = function(x)
{
   p = ncol(x)
   n = nrow(x)
   center = colMeans(x)
   x.mean = x - matrix(rep(center, n), n, p, byrow=T)
   scale = sqrt(colSums(x.mean^2)/n)
   xx = t(t(x.mean)/scale)
   list(xx, center, scale)
}

