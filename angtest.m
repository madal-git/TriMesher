function ang = angle(u,v)
nu = norm(u);
nv = norm(v);

numer = u'*v;
costheta = numer/(nu*nv);

ang = acosd(costheta);
