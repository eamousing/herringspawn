bilin_inv <- function(f, g, Farr, Garr, imax, jmax)
{
  b <- vector(length = 2)
  maxiter <- 10
  
  b[1] <- imax/2
  b[2] <- jmax/2
  for(iter in 1:maxiter)
  {
    i <- as.integer(b[1])
    j <- as.integer(b[2])
    i <- max(2, min(i, imax - 1))
    j <- max(2, min(j, jmax - 1))
    p <- b[1] - i
    q <- b[2] - j
    
    fs <- (1-p) * (1-q) * Farr[i,j] + p*(1-q) * Farr[i+1, j] +
      (1-p) * q * Farr[i,j+1] + p*q*Farr[i+1,j+1]
    gs <- (1-p) * (1-q) * Garr[i,j] + p*(1-q) * Garr[i+1, j] +
      (1-p) * q * Garr[i,j+1] + p*q*Garr[i+1,j+1]
    if(((fs-f)^2 + (gs-g)^2) > 1e-12)
    {
      Fx <- (1-q)*(Farr[i+1,j]-Farr[i,j]) + q*(Farr[i+1,j+1]-Farr[i,j+1])
      Fy <- (1-p)*(Farr[i,j+1]-Farr[i,j]) + p*(Farr[i+1,j+1]-Farr[i+1,j])
      Gx <- (1-q)*(Garr[i+1,j]-Garr[i,j]) + q*(Garr[i+1,j+1]-Garr[i,j+1])
      Gy <- (1-p)*(Garr[i,j+1]-Garr[i,j]) + p*(Garr[i+1,j+1]-Garr[i+1,j])
      
      det <- Fx*Gy - Fy*Gx
      
      if(det != 0)
      {
        b[1] = b[1] - (Gy*(fs-f) - Fy*(gs-g)) / det
        b[2] = b[2] - (-Gx*(fs-f) + Fx*(gs-g)) / det
      }
    }
    else
    {
      if(i >= imax-1 | i <= 2 | j >= jmax-1 | j <= 2 & 
         b[1] >= imax-1 | b[1] <= 0 | b[2] >= jmax-1 | b[2] <= 2)
      {
        return(NA)
      }
      return(b)
    }
  }
}