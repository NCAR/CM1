NCAR_primary-stylelinks arc_iframe-stylelinks:
	topdir=$$(pwd) ; \
	target=$(subst -stylelinks,,$@) ; \
	cd $${topdir} && ln -sf themes/$${target}.yml theme.yml
