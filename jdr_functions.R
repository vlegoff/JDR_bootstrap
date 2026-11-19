prep_data = function(data, instance=NULL) {
    if (!is.null(instance)) {
        data = data[instance,]
    }

    blocks = c("mutation", "rna", "mirna", "cnv", "clinical")
    x_list = lapply(blocks,
        \(bl) data[,.SD, .SDcols=grepl(paste0("_", bl), colnames(data))])
    names(x_list) = blocks
    x_list = lapply(x_list, as.matrix)
    clin_data = x_list[["clinical"]]
    x_list = x_list[c("mutation", "rna", "mirna", "cnv")]

    time = data$time
    status = data$status
    ysurv = survival::Surv(time=time, event=status)
    return(list(x=x_list, clin=clin_data, y=ysurv))
}

rgcca_wrapper = function(data, job, instance, clin, sup, n, tau, lambda,
    weights_orig=NULL) {
    #source("jdr_functions.R") # for batchtools...
    #data = prep_data(data, instance)
    x = instance$x
    blocks = names(x)
    clin_data = instance$clin
    y = instance$y
    weights_orig = weights_orig$jdr
    
    response = NULL
    if(sup) {
        null_mod = survival::coxph(y~1)
        x[["residuals"]] = residuals(null_mod, type="deviance")
        response = length(x)
    } else {
        complete_matrix = matrix(1, length(x),
                               length(x))
        diag(complete_matrix) = 0
        rownames(complete_matrix) = names(x)
        colnames(complete_matrix) = names(x)
    }
    rgcca_args = list(
        response=response,
        connection=if(!sup) complete_matrix,
        tau=tau,
        ncomp=n,
        scheme="factorial",
        init="svd", 
        method="rgcca",
        scale=TRUE,
        scale_block="inertia",
        verbose=F
    )
    rgcca_fit = mlr3misc::invoke(RGCCA::rgcca,
        .args=c(list(blocks=x), rgcca_args))
    comps = Reduce(cbind,
        rgcca_fit$Y[names(rgcca_fit$Y)!="residuals"])
    coln = paste(
        rep(blocks, each=rgcca_args$ncomp),
        rep(1:rgcca_args$ncomp, length(blocks)),
        sep="."
    )
    colnames(comps) = coln

    weights = rgcca_fit$a
    if (!is.null(weights_orig)) {
        vars_bts = lapply(weights, rownames)
        vars_orig = lapply(weights_orig, rownames)
        common_vars = mapply(intersect, vars_bts, vars_orig, SIMPLIFY=FALSE)
        common_vars = common_vars[names(common_vars)!="residuals"] # to remove residuals
        for (b in blocks) {
            weights[[b]] = as.matrix(weights[[b]][common_vars[[b]],])
            weights_orig[[b]] = as.matrix(weights_orig[[b]][common_vars[[b]],])
        }
        cors = lapply(blocks, \(b) diag(cor(weights[[b]],
            weights_orig[[b]])))
        names(cors) = blocks
        cors = lapply(cors, sign)
        print(cors)
        #weights = lapply(blocks, \(b) weights[[b]]*cors[[b]]) # wrong
        #weights = lapply(blocks, \(b) sweep(weights[[b]], 2, cors[[b]], "*"))
        weights = lapply(blocks, \(b) t(t(weights[[b]])*cors[[b]]))
        names(weights) = blocks
        comps_old = comps
        #comps = comps * Reduce(c, cors)
        comps = sweep(comps, 2, Reduce(c, cors), "*")
        #print(Reduce(c, cors))
        print(diag(cor(comps, comps_old)))
    }

    if(clin) {
        comps = cbind(comps, clin_data)
    }

    glmnet_fit = glmnet::glmnet(comps, y, family="cox", alpha=0,
        lambda=lambda)
    sigma = apply(comps, 2, sd)
    beta = glmnet_fit$beta[,1]
    beta = beta*sigma

    learner_id = paste0("rgcca_", ifelse(sup, "sup", "unsup"),
        "_", ifelse(clin, "clin_fav", ""))
    return(list(jdr=weights, cox=beta, model="rgcca", learner_id=learner_id,
                latent_space=comps))
}

jive_wrapper = function(data, job, instance, clin, j, a, lambda, weights_orig=NULL){
    #data = prep_data(data, instance)
    #source("jdr_functions.R")
    x = instance$x
    blocks = names(x)
    clin_data = instance$clin
    y = instance$y
    weights_orig = weights_orig$jdr

    # normalization
    n = sapply(x, function(xi) prod(dim(xi)))
    x = lapply(x, scale, center=TRUE, scale=FALSE)
    x = lapply(x, \(xi) xi/(norm(xi, "F")*sqrt(sum(n))))
    x = lapply(x, t) # why oh god why

    if (a==0) {
        x = Reduce(rbind, x)
        jive_fit = svd(x, nu=j, nv=j)
        joint.load = jive_fit$u
        if (!is.null(weights_orig)) {
            svd_tmp = svd(t(weights_orig$joint)%*%joint.load)
            R = svd_tmp$u %*% t(svd_tmp$v)
            joint.load = joint.load %*% t(R)
        }
        weights = list(joint=joint.load, indiv=NULL)
        rownames(weights$joint) = rownames(x)
        latent_space = t(x) %*% joint.load
        #latent_space = jive_fit$v
        #print(dim(latent_space))
        #if (!is.null(weights_orig)) latent_space = latent_space %*% t(R)
        colnames(latent_space) = paste0("common.", 1:j)
    }
    else if (j==0) {
        jive_fit = lapply(x, svd, nu=a, nv=a)
        indiv.load = lapply(jive_fit, "[[", "u")
        if (!is.null(weights_orig)) {
            for (i in seq_along(blocks)) {
                svd_tmp = svd(t(weights_orig$indiv[[i]])%*%indiv.load[[i]]) 
                R = svd_tmp$u %*% t(svd_tmp$v)
                indiv.load[[i]] = indiv.load[[i]] %*% t(R)
            }
            names(indiv.load) = blocks
        }
        weights = list(joint=NULL, indiv=indiv.load)
        for (i in seq_along(blocks)) {
            rownames(weights$indiv[[i]]) = rownames(x[[i]])
        }
        names(weights$indiv) = blocks
        #latent_space = Reduce(cbind, mapply(\(xi, wi) t(xi) %*% wi, x, indiv.load))
        #latent_space = .mapply(\(x, w) {t(x)%*%w}, dots=list(x=x, w=indiv.load))
        latent_space = mapply(\(x, w) {t(x)%*%w}, x, indiv.load, SIMPLIFY=FALSE)
        latent_space = Reduce(cbind, latent_space)
        colnames(latent_space) = paste(
            rep(names(x), each=a),
            rep(1:a, length(names(x))),
            sep="."
        )
    } else {
        jive_args = list(
            rankJ=j,
            rankA=rep(a, length(x)),
            method="given",
            orthIndiv=FALSE,
            scale=FALSE,
            center=FALSE
        )
        jive_fit = mlr3misc::invoke(r.jive::jive,
            .args=c(list(data=x), jive_args))
        predict = r.jive::jive.predict(data.new=x, jive_fit)

        # rotation
        weights = list(joint=predict$joint.load, indiv=predict$indiv.load)
        joint.scores = t(predict$joint.scores)
        indiv.scores = lapply(predict$indiv.scores, t)
        if (!is.null(weights_orig)){
            # joint loadings
            svd_tmp = svd(t(weights_orig$joint)%*%weights$joint)
            R = svd_tmp$u %*% t(svd_tmp$v)
            weights$joint = weights$joint %*% t(R)
            joint.scores = t(predict$joint.scores) %*% t(R)
            # indiv loadings
            for (i in seq_along(blocks)) {
                svd_tmp = svd(t(weights_orig$indiv[[i]])%*%weights$indiv[[i]]) 
                R = svd_tmp$u %*% t(svd_tmp$v)
                weights$indiv[[i]] = weights$indiv[[i]] %*% t(R)
                indiv.scores[[i]] = t(predict$indiv.scores[[i]]) %*% t(R)
            }
            names(weights$indiv) = blocks
        }
        rownames(weights$joint) = Reduce(c, lapply(x, rownames))
        for (i in seq_along(blocks)) {
            rownames(weights$indiv[[i]]) = rownames(x[[i]])
        }
        names(weights$indiv) = blocks

        latent_space = matrix(0,
            nrow=ncol(x[[1]]),
            ncol=jive_args$rankJ+length(x)*jive_args$rankA)

        latent_space[,1:jive_args$rankJ] = joint.scores
        
        latent_space[,(jive_args$rankJ+1):ncol(latent_space)] = 
            Reduce(cbind, indiv.scores)

        colnames(latent_space) = c(
            paste0("common.", 1:j),
            paste(
                rep(names(x), each=a),
                rep(1:a, length(names(x))),
                sep="."
            )
        )

    }

    if (clin) {
        latent_space = cbind(latent_space, clin_data)
    }

    if (ncol(latent_space)==1) {
        cox_fit = survival::coxph(y ~ latent_space, ties="breslow")
        beta = coef(cox_fit)
        names(beta) = names(latent_space) # to get common.1 as name
    } else {
        glmnet_fit = glmnet::glmnet(latent_space, y, family="cox", alpha=0,
            lambda=lambda)
        beta = glmnet_fit$beta[,1]
    }
    sigma = apply(latent_space, 2, sd)
    beta = beta*sigma

    learner_id = paste0("jive_", ifelse(clin, "clin_fav", ""))
    return(list(jdr=weights, cox=beta, model="jive", learner_id=learner_id,
                latent_space=latent_space))
}

ipca_wrapper = function(data, job, instance, clin, a, lambda, weights_orig=NULL){
    #data = prep_data(data, instance)
    #source("jdr_functions.R")
    x = instance$x
    blocks = names(x)
    clin_data = instance$clin
    y = instance$y
    weights_orig = weights_orig$jdr

    # normalization
    cols = lapply(x, function(xi) apply(xi, 2, stats::var)!=0)
    x = lapply(names(x), function(n) x[[n]][,cols[[n]]])
    x = lapply(x, scale, center=TRUE, scale=TRUE)
    #x = lapply(x, \(xi) xi/norm(xi, "F"))

    pca_fit = lapply(x, svd, nu=a, nv=a)
    weights = lapply(pca_fit, "[[", "v")
    for (i in seq_along(blocks)) rownames(weights[[i]]) = colnames(x[[i]])
    if (!is.null(weights_orig)) {
        vars_bts = lapply(weights, rownames)
        vars_orig = lapply(weights_orig, rownames)
        common_vars = mapply(intersect, vars_bts, vars_orig, SIMPLIFY=FALSE)
        names(common_vars) = blocks
        for (i in seq_along(blocks)) {
            weights_tmp = as.matrix(weights[[i]][common_vars[[i]],])
            weights_orig_tmp = as.matrix(weights_orig[[i]][common_vars[[i]],])
            svd_tmp = svd(t(weights_orig_tmp)%*%weights_tmp) 
            R = svd_tmp$u %*% t(svd_tmp$v)
            weights[[i]] = weights[[i]] %*% t(R)
        }
        names(weights) = blocks
    }
    for (i in seq_along(blocks)) {
        rownames(weights[[i]]) = colnames(x[[i]])
    }
    names(weights) = blocks
    latent_space = mapply(\(x, w) {x%*%w}, x, weights, SIMPLIFY=FALSE)
    latent_space = Reduce(cbind, latent_space)
    colnames(latent_space) = paste(
        rep(blocks, each=a),
        rep(1:a, length(blocks)),
        sep="."
    )

    if (clin) {
        latent_space = cbind(latent_space, clin_data)
    }

    if (ncol(latent_space)==1) {
        cox_fit = survival::coxph(y ~ latent_space, ties="breslow")
        beta = coef(cox_fit)
        names(beta) = names(latent_space) # to get common.1 as name
    } else {
        glmnet_fit = glmnet::glmnet(latent_space, y, family="cox", alpha=0,
            lambda=lambda)
        beta = glmnet_fit$beta[,1]
    }
    sigma = apply(latent_space, 2, sd)
    beta = beta*sigma

    learner_id = paste0("ipca_", ifelse(clin, "clin_fav", ""))
    return(list(jdr=weights, cox=beta, model="ipca", learner_id=learner_id,
                latent_space=latent_space))
}

sjive_wrapper = function(data, job, instance, clin, j, a, eta, lambda, weights_orig=NULL){
    #data = prep_data(data, instance)
    #source("jdr_functions.R")
    x = instance$x
    blocks = names(x)
    clin_data = instance$clin
    y = instance$y
    weights_orig = weights_orig$jdr

    # normalization
    cols = lapply(x, function(xi) apply(xi, 2, stats::var)!=0)
    x = lapply(names(x), function(n) x[[n]][,cols[[n]]])
    x = lapply(x, scale, center=TRUE, scale=TRUE)
    #x = lapply(x, \(xi) xi/norm(xi, "F"))
    x = lapply(x, t) # why oh god why

    null_mod = survival::coxph(y~1)
    resi = residuals(null_mod, type="deviance")
    resi = as.numeric(scale(as.numeric(resi)))

    jive_args = list(
        rankJ=j,
        rankA=rep(a, length(x)),
        eta=eta,
        center.scale=FALSE,
        reduce.dim=TRUE
    )
    jive_fit = mlr3misc::invoke(sup.r.jive::sJIVE,
        .args=c(list(X=x, Y=resi), jive_args))
    predict = sup.r.jive:::predict.sJIVE(jive_fit, newdata=x)

    # rotation
    weights = list(joint=Reduce(rbind, jive_fit$U_I), indiv=jive_fit$W_I)
    rownames(weights$joint) = unlist(sapply(x, rownames))
    weights$indiv = mapply(
        function(w, xi) {
            rownames(w) = rownames(xi)
            return(w)
        },
        weights$indiv, x)
    joint.scores = t(predict$Sj)
    indiv.scores = lapply(predict$Si, t)
    if (!is.null(weights_orig)){
        message("Aligning...")
        print(names(weights_orig))
        print(sapply(weights_orig, is.null))
        vars_bts = lapply(weights$indiv, rownames)
        vars_orig = lapply(weights_orig$indiv, rownames)
        common_vars = mapply(intersect, vars_bts, vars_orig, SIMPLIFY=FALSE)
        names(common_vars) = blocks
        common_vars_joint = Reduce(c, common_vars)
        # joint loadings
        weights_tmp = as.matrix(weights$joint[common_vars_joint,])
        weights_orig_tmp = as.matrix(weights$joint[common_vars_joint,])
        svd_tmp = svd(t(weights_orig_tmp)%*%weights_tmp) 
        R = svd_tmp$u %*% t(svd_tmp$v)
        weights$joint = weights$joint %*% t(R)
        joint.scores = t(predict$Sj) %*% t(R)
        # indiv loadings
        for (i in seq_along(blocks)) {
            weights_tmp = as.matrix(weights$indiv[[i]][common_vars[[i]],])
            weights_orig_tmp = as.matrix(weights_orig$indiv[[i]][common_vars[[i]],])
            svd_tmp = svd(t(weights_orig_tmp)%*%weights_tmp) 
            R = svd_tmp$u %*% t(svd_tmp$v)
            weights$indiv[[i]] = weights$indiv[[i]] %*% t(R)
            indiv.scores[[i]] = t(predict$Si[[i]]) %*% t(R)
        }
        names(weights$indiv) = blocks
    }
    rownames(weights$joint) = Reduce(c, lapply(x, rownames))
    for (i in seq_along(blocks)) {
        rownames(weights$indiv[[i]]) = rownames(x[[i]])
    }
    names(weights$indiv) = blocks
    names(x) = blocks

    latent_space = matrix(
        0,
        nrow=ncol(x[[1]]),
        ncol=jive_args$rankJ+length(x)*jive_args$rankA[1]
    )
    if(jive_args$rankA[1]==0) {
        latent_space[,1:jive_args$rankJ] = t(joint.scores)
        colnames(latent_space) = paste0("common.", 1:jive_args$rankJ)
    } else if (jive_args$rankJ==0) {
        latent_space[,1:ncol(latent_space)] = 
            t(Reduce(rbind, indiv.scores)) # cbind was not working

        colnames(latent_space) =
            paste0(
                rep(names(x), each=jive_args$rankA[1]),
                ".",
                rep(1:jive_args$rankA[1], length(names(x)))
            )
    } else {
        latent_space[,1:jive_args$rankJ] = t(joint.scores)
      
        latent_space[,(jive_args$rankJ+1):ncol(latent_space)] = 
            t(Reduce(rbind, indiv.scores)) # cbind was not working
        colnames(latent_space) = c(
            paste0("common.", 1:jive_args$rankJ),
            paste0(
                rep(names(x), each=jive_args$rankA[1]),
                ".",
                rep(1:jive_args$rankA[1], length(names(x)))
            )
        )
    }

    if (clin) {
        latent_space = cbind(latent_space, clin_data)
    }

    if (ncol(latent_space)==1) {
        cox_fit = survival::coxph(y ~ latent_space, ties="breslow")
        beta = coef(cox_fit)
        names(beta) = names(latent_space) # to get common.1 as name
    } else {
        glmnet_fit = glmnet::glmnet(latent_space, y, family="cox", alpha=0,
            lambda=lambda)
        beta = glmnet_fit$beta[,1]
    }
    sigma = apply(latent_space, 2, sd)
    beta = beta*sigma

    learner_id = paste0("sjive_", ifelse(clin, "clin_fav", ""))
    return(list(jdr=weights, cox=beta, model="sjive", learner_id=learner_id))
}
