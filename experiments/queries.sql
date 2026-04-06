SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.mol FROM compound_structures m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE is_substruct(m.mol, 'COc1cc2c(Nc3cc(CC(=O)Nc4cccc(F)c4F)[nH]n3)ncnc2cc1OCCCN(CCO)CC(C)C');

SELECT count(*) FROM compound_structures m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      WHERE is_substruct(m.mol, 'O=CNCCc1ccccc1');

SELECT avg(a.value), stddev(a.value), a.units,a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.mol FROM compound_structures m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE is_substruct(m.mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
      GROUP BY m.mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.mol FROM compound_structures m
    INNER JOIN activities a ON a.molregno=m.molregno
    INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
    INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
    INNER JOIN assays ys ON ys.tid=bs.tid
    WHERE is_substruct(m.mol, 'N1C=CC=N1')
    GROUP BY m.mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

