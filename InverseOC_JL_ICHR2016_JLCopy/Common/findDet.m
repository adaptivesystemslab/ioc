function determ = findDet(H)
  e1 = H(1,1);
  e2 = H(1,2);
  e3 = H(1,3);
  
  s1 = e1 * (H(2,2)*H(3,3) - H(3,2)*H(2,3));
  s2 = -e2 * (H(2,1)*H(3,3) - H(3,1)*H(2,3));
  s3 = e3 * (H(2,1)*H(3,2) - H(3,1)*H(2,2));
  
  y1 = simplify(s1);
  y2 = simplify(s2);
  y3 = simplify(s3);
  
  determ = s1 + s2 + s3;
end