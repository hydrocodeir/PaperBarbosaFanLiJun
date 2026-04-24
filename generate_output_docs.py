from __future__ import annotations

import os
import re
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"
FIGURES = OUTPUTS / "figures"
TABLES = OUTPUTS / "tables"
DOCS = OUTPUTS / "output_docs"


def tau_label(value: float | str) -> str:
    return f"{float(value):0.2f}"


def index_title_fa(index_name: str) -> str:
    mapping = {
        "warm_days": "روزهای گرم",
        "warm_nights": "شب‌های گرم",
        "cool_days": "روزهای سرد",
        "cool_nights": "شب‌های سرد",
    }
    return mapping.get(index_name, index_name.replace("_", " "))


def index_title_en(index_name: str) -> str:
    mapping = {
        "warm_days": "Warm Days",
        "warm_nights": "Warm Nights",
        "cool_days": "Cool Days",
        "cool_nights": "Cool Nights",
    }
    return mapping.get(index_name, index_name.replace("_", " ").title())


def make_rel(from_file: Path, target: Path) -> str:
    return os.path.relpath(target, start=from_file.parent).replace("\\", "/")


def write_doc(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.strip() + "\n", encoding="utf-8")


def build_context(cfg: dict) -> dict:
    start_year, end_year = cfg["data"]["analysis_years"]
    n_years = end_year - start_year + 1
    first_end = start_year + (n_years // 2) - 1
    second_start = first_end + 1
    sample_paper2_fig1 = sorted((FIGURES / "paper2_station_figures" / "figure1_timeseries").glob("*.png"))
    sample_paper2_fig2 = sorted((FIGURES / "paper2_station_figures" / "figure2_quantile_coefficients").glob("*.png"))
    sample_paper1_fig4 = sorted((FIGURES / "paper1_station_figures" / "figure4_bootstrap_distributions").glob("*.png"))
    return {
        "analysis_years": f"{start_year}-{end_year}",
        "focus_quantiles": [float(x) for x in cfg["quantile_regression"]["focus_quantiles"]],
        "focus_quantiles_text": ", ".join(tau_label(x) for x in cfg["quantile_regression"]["focus_quantiles"]),
        "full_quantile_start": tau_label(cfg["quantile_regression"]["full_quantiles"]["start"]),
        "full_quantile_stop": tau_label(cfg["quantile_regression"]["full_quantiles"]["stop"]),
        "full_quantile_step": cfg["quantile_regression"]["full_quantiles"]["step"],
        "bootstrap_method": cfg["bootstrap"]["method"],
        "bootstrap_reps": cfg["bootstrap"]["n_reps"],
        "bootstrap_deeper_reps": cfg["advanced_analyses"]["bootstrap_depth_sensitivity"]["n_reps"],
        "primary_delta": cfg["feature_engineering"]["primary_delta"],
        "cluster_algorithm": cfg["clustering"]["algorithm"],
        "cluster_linkage": cfg["clustering"]["linkage"],
        "cluster_metric": cfg["clustering"]["metric"],
        "cluster_count": cfg["clustering"]["n_clusters"],
        "interp_method": cfg["spatial_visualization"]["interpolation_method"],
        "interp_smooth": cfg["spatial_visualization"]["interpolation_smooth"],
        "ref_low": cfg["index_construction"]["lower_percentile"],
        "ref_high": cfg["index_construction"]["upper_percentile"],
        "window_days": cfg["index_construction"]["percentile_window_days"],
        "min_coverage": cfg["index_construction"]["annual_min_valid_coverage_pct"],
        "split_1": f"{start_year}-{first_end}",
        "split_2": f"{second_start}-{end_year}",
        "paper2_fig1_sample": sample_paper2_fig1[0] if sample_paper2_fig1 else None,
        "paper2_fig2_sample": sample_paper2_fig2[0] if sample_paper2_fig2 else None,
        "paper1_fig4_sample": sample_paper1_fig4[0] if sample_paper1_fig4 else None,
    }


def bilingual_paragraph(fa: str, en: str) -> str:
    return "\n".join(
        [
            "## فارسی",
            "",
            fa,
            "",
            "## English",
            "",
            en,
        ]
    )


def figure_doc(path: Path, figure_path: Path, desc: dict[str, str]) -> None:
    rel = make_rel(path, figure_path)
    content = f"""
# `{figure_path.name}`

**نوع / Type:** {desc["kind_fa"]} / {desc["kind_en"]}

[لینک مستقیم به شکل / Direct figure link]({rel})

![{figure_path.name}]({rel})

{bilingual_paragraph(desc["fa"], desc["en"])}
"""
    write_doc(path, content)


def table_doc(path: Path, table_path: Path, desc: dict[str, str]) -> None:
    rel = make_rel(path, table_path)
    content = f"""
# `{table_path.name}`

[لینک مستقیم به جدول / Direct table link]({rel})

{bilingual_paragraph(desc["fa"], desc["en"])}
"""
    write_doc(path, content)


def fallback_table_desc(name: str) -> dict[str, str]:
    return {
        "fa": f"این جدول یکی از خروجی‌های ساختاریافته‌ی پایپ‌لاین است و اطلاعات مرحله‌ای مرتبط با فایل `{name}` را برای مستندسازی بازتولیدپذیری و رهگیری فرآیند تحلیل نگه می‌دارد.",
        "en": f"This table is a structured pipeline export that preserves the stage-specific information associated with `{name}` for reproducibility and workflow traceability.",
    }


def describe_table(table_path: Path, ctx: dict) -> dict[str, str]:
    name = table_path.name
    p = ctx["primary_delta"]
    fq = ctx["focus_quantiles_text"]
    descriptions: dict[str, dict[str, str]] = {
        "data_quality_station_summary.csv": {
            "fa": "این جدول خلاصه‌ی کنترل کیفیت در مقیاس ایستگاه را ارائه می‌کند و برای هر ایستگاه درصد کامل‌بودن `Tmin`، `Tmax` و `Tmean`، تعداد تاریخ‌های تکراری، و ناسازگاری‌های اولیه‌ی دمایی را ثبت می‌نماید. این فایل در مرحله‌ی ممیزی داده تولید شده و نقش آن معرفی رسمی کیفیت داده‌های ورودی پیش از ساخت شاخص‌های سالانه است.",
            "en": "This table provides the station-level data-quality summary, including completeness percentages for `Tmin`, `Tmax`, and `Tmean`, counts of duplicate dates, and basic temperature-consistency checks. It is produced during the data-audit stage and serves as a formal description of input-data quality prior to annual-index construction.",
        },
        "data_homogeneity_tests_station_summary.csv": {
            "fa": "این جدول نتایج آزمون‌های همگنی را در سطح ایستگاه نگه می‌دارد و آماره‌ها، `p-value`ها و سال تغییر محتمل را برای آزمون‌های `Pettitt`، `SNHT` و `Buishand` در نسخه‌های خام و روندزدایی‌شده گزارش می‌کند. هدف آن مستندسازی کامل ارزیابی شکست ساختاری و ناهمگنی احتمالی در سری‌های سالانه است.",
            "en": "This table stores the station-level homogeneity-test results, including statistics, `p-values`, and candidate change years for the `Pettitt`, `SNHT`, and `Buishand` tests under both raw and detrended formulations. Its purpose is to document the structural-break and homogeneity screening applied to the annual series.",
        },
        "data_quality_homogeneity_overview.csv": {
            "fa": "این جدول یک نمای فشرده از کیفیت و همگنی شبکه فراهم می‌کند و شاخص‌های خلاصه‌ای مانند تعداد ایستگاه‌ها، میانه‌ی کامل‌بودن داده‌ها، تعداد ایستگاه‌های دارای داده‌ی تکراری یا ناسازگار، و شمار ایستگاه‌های پرچم‌خورده در آزمون‌های همگنی را یکجا جمع می‌کند. این فایل برای ارائه‌ی یک مرور رسمی و سریع از وضعیت داده‌های ورودی مناسب است.",
            "en": "This table provides a compact network-level overview of data quality and homogeneity, bringing together summary indicators such as the number of stations, median completeness, counts of duplicated or inconsistent records, and numbers of stations flagged by homogeneity tests. It is intended as a formal high-level synopsis of the input-data status.",
        },
        "annual_extreme_indices.csv": {
            "fa": f"این جدول خروجی اصلی مرحله‌ی ساخت شاخص‌های سالانه است و برای هر ایستگاه و هر سال، شمار سالانه‌ی شاخص‌های `warm_days`، `warm_nights`، `cool_days` و `cool_nights` را ثبت می‌کند. این شاخص‌ها با استفاده از آستانه‌های صدکی `{ctx['ref_high']}` و `{ctx['ref_low']}`، پنجره‌ی تقویمی `{ctx['window_days']}` روزه، و شرط حداقل `{ctx['min_coverage']}` درصد پوشش معتبر سالانه تولید شده‌اند.",
            "en": f"This table is the primary output of annual-index construction and records yearly counts of `warm_days`, `warm_nights`, `cool_days`, and `cool_nights` for each station and year. These indices were derived using the `{ctx['ref_high']}`th and `{ctx['ref_low']}`th percentile thresholds, a `{ctx['window_days']}`-day moving calendar window, and a minimum annual valid-coverage requirement of `{ctx['min_coverage']}` percent.",
        },
        "qr_all_quantiles_long.csv": {
            "fa": f"این جدول خروجی کامل رگرسیون صدکی بر روی شبکه‌ی متراکم صدک‌ها است و برآوردهای شیب را برای کل بازه‌ی `τ` از `{ctx['full_quantile_start']}` تا `{ctx['full_quantile_stop']}` با گام `{ctx['full_quantile_step']}` نگه می‌دارد. این فایل مرجع اصلی برای ترسیم پروفایل‌های صدکی پیوسته و شکل‌هایی است که رفتار روند را در سراسر توزیع نمایش می‌دهند.",
            "en": f"This table contains the full quantile-regression output across the dense quantile grid, preserving slope estimates for the complete `τ` range from `{ctx['full_quantile_start']}` to `{ctx['full_quantile_stop']}` with step `{ctx['full_quantile_step']}`. It serves as the main source for continuous quantile-profile graphics and for figures that display trend behavior across the distribution.",
        },
        "qr_focus_slopes_and_bootstrap_summary.csv": {
            "fa": f"این جدول خلاصه‌ی اصلی رگرسیون صدکی در صدک‌های کانونی `{fq}` است و برای هر ایستگاه و هر شاخص، شیب‌های صدکی، بازه‌های اطمینان تحلیلی، و خلاصه‌های بوت‌استرپ را گردآوری می‌کند. این فایل مبنای اصلی شکل‌های ایستگاهی، نقشه‌های شیب، و اسناد عدم‌قطعیت در مستندات تکمیلی است.",
            "en": f"This table is the principal quantile-regression summary at the focal quantiles `{fq}`, compiling quantile slopes, analytic confidence intervals, and bootstrap summaries for each station and index. It provides the core input for station-level graphics, slope maps, and uncertainty-oriented supplementary documentation.",
        },
        "clustering_feature_table.csv": {
            "fa": f"این جدول ماتریس ویژگی مورد استفاده در مرحله‌ی خوشه‌بندی را نگه می‌دارد و شیب‌های صدکی، شاخص‌های مشتق‌شده مانند `{p}`، خلاصه‌های بوت‌استرپ، و برچسب‌های خوشه را در یک قاب تحلیلی واحد ترکیب می‌کند. از این رو، این فایل نمایش عددی امضای روند هر ایستگاه برای منطقه‌بندی و نمایش‌های مکمل فضایی و مقایسه‌ای است.",
            "en": f"This table stores the feature matrix used in the clustering stage, combining quantile slopes, derived metrics such as `{p}`, bootstrap summaries, and cluster labels within a single analytical frame. It therefore represents the numerical trend signature of each station for regionalization and related spatial or comparative displays.",
        },
        "clustering_feature_screening_summary.csv": {
            "fa": "این جدول فرایند غربال‌گری ویژگی‌ها را پیش از خوشه‌بندی مستند می‌کند و نشان می‌دهد کدام متغیرها به‌دلیل هم‌بستگی بسیار بالا حذف یا حفظ شده‌اند. نقش آن شفاف‌سازی این است که مجموعه‌ی نهایی ویژگی‌ها چگونه برای پرهیز از افزونگی آماری انتخاب شده است.",
            "en": "This table documents the pre-clustering feature-screening procedure by showing which variables were retained or removed because of very high pairwise correlation. Its role is to clarify how the final feature set was selected to reduce statistical redundancy.",
        },
        "cluster_assignments.csv": {
            "fa": f"این جدول برچسب خوشه‌ی پایه را برای هر ایستگاه و هر شاخص ارائه می‌کند. این تخصیص‌ها با خوشه‌بندی `{ctx['cluster_algorithm']}` تحت تنظیمات `linkage={ctx['cluster_linkage']}`، `metric={ctx['cluster_metric']}` و `n_clusters={ctx['cluster_count']}` به‌دست آمده‌اند و مبنای نقشه‌های خوشه و انتخاب ایستگاه‌های نماینده را تشکیل می‌دهند.",
            "en": f"This table reports the baseline cluster label for each station and index. These assignments were obtained using `{ctx['cluster_algorithm']}` clustering with `linkage={ctx['cluster_linkage']}`, `metric={ctx['cluster_metric']}`, and `n_clusters={ctx['cluster_count']}`, and they form the basis for cluster maps and representative-station selection.",
        },
        "cluster_assignments_reduced_features.csv": {
            "fa": "این جدول خروجی اجرای حساسیت خوشه‌بندی با مجموعه‌ی دیگری از ویژگی‌ها را ارائه می‌کند. کارکرد آن این است که امکان مقایسه‌ی مستقیم بین منطقه‌بندی پایه و اجرای حساسیت فراهم شود و پایداری ساختار خوشه‌ها به‌صورت شفاف مستند گردد.",
            "en": "This table reports the clustering-sensitivity rerun based on an alternative feature set. Its function is to enable direct comparison between the baseline regionalization and the sensitivity configuration, thereby documenting cluster stability explicitly.",
        },
        "cluster_robustness_summary.csv": {
            "fa": "این جدول خلاصه‌ی پایداری خوشه‌بندی را با معیارهایی نظیر `Adjusted Rand Index` و سهم برچسب‌های یکسان ارائه می‌کند. بدین‌ترتیب، فایل حاضر برای توصیف رسمی میزان سازگاری ساختار خوشه‌ها بین اجراهای مختلف خوشه‌بندی به‌کار می‌رود.",
            "en": "This table summarizes clustering robustness using metrics such as the `Adjusted Rand Index` and the fraction of matching labels. It is therefore used to formally describe the degree of agreement in cluster structure across different clustering runs.",
        },
        "alternative_clustering_sensitivity_summary.csv": {
            "fa": "این جدول مقایسه‌ی روش پایه‌ی خوشه‌بندی با روش‌های جایگزین را خلاصه می‌کند و میزان توافق ساختار خوشه‌ها را برای الگوریتم‌ها و تنظیمات مختلف گزارش می‌نماید. هدف آن مستندسازی حساسیت منطقه‌بندی نسبت به انتخاب روش خوشه‌بندی است.",
            "en": "This table summarizes the comparison between the baseline clustering strategy and alternative methods, reporting the degree of agreement in cluster structure across multiple algorithms and configurations. Its purpose is to document the sensitivity of regionalization to clustering choice.",
        },
        "alternative_clustering_assignments.csv": {
            "fa": "این جدول برچسب خوشه‌ی هر ایستگاه را برای روش پایه و همه‌ی روش‌های جایگزین در کنار یکدیگر نگه می‌دارد. این چیدمان امکان رهگیری تغییرات جایگاه ایستگاه‌ها را در سناریوهای مختلف خوشه‌بندی فراهم می‌سازد.",
            "en": "This table stores the cluster label of each station under the baseline and all alternative clustering methods side by side. Such a layout enables direct tracking of how station membership changes across clustering scenarios.",
        },
        "publication_summary_table.csv": {
            "fa": "این جدول نسخه‌ی فشرده و ارائه‌محور خروجی‌های اصلی است و فقط متغیرهای کلیدی مورد نیاز برای گزارش‌نویسی، شکل‌ها و جمع‌بندی نتایج را نگه می‌دارد. از این رو، این فایل را می‌توان جدول مرجع سبک‌وزن برای استفاده در متن مقاله و پیوست‌ها دانست.",
            "en": "This table is a compact, presentation-oriented summary of the main outputs and retains only the key variables needed for reporting, figure preparation, and result synthesis. It may therefore be regarded as a lightweight reference table for use in the manuscript and appendices.",
        },
        "homogeneity_flag_exclusion_sensitivity.csv": {
            "fa": "این جدول حساسیت خلاصه‌های منطقه‌ای را نسبت به حذف ایستگاه‌های پرچم‌خورده در آزمون همگنی مستند می‌کند. در این فایل، اختلاف بین اجرای کامل و اجرای بدون ایستگاه‌های مسئله‌دار برای شیب‌ها و معیارهای مشتق‌شده ثبت شده است.",
            "en": "This table documents the sensitivity of regional summaries to the exclusion of stations flagged by homogeneity screening. It records the differences between the full analysis and the reduced analysis for slopes and derived metrics.",
        },
        "bootstrap_distributions_long.csv": {
            "fa": f"این جدول توزیع کامل بوت‌استرپ را در قالب بلند نگه می‌دارد و برای هر تکرار، ایستگاه، شاخص و صدک، برآوردهای شیب را ثبت می‌کند. بوت‌استرپ اصلی با روش `{ctx['bootstrap_method']}` و `{ctx['bootstrap_reps']}` تکرار اجرا شده است و این فایل ماده‌ی خام لازم برای شکل‌های توزیع بوت‌استرپ و توصیف عدم‌قطعیت را فراهم می‌سازد.",
            "en": f"This table preserves the full bootstrap distribution in long format, recording slope estimates for each replicate, station, index, and quantile. The main bootstrap was run with the `{ctx['bootstrap_method']}` method and `{ctx['bootstrap_reps']}` replicates, making this file the raw source for bootstrap-distribution figures and uncertainty documentation.",
        },
        "bootstrap_depth_sensitivity_summary.csv": {
            "fa": f"این جدول مقایسه‌ی خلاصه‌ای بین بوت‌استرپ پایه با `{ctx['bootstrap_reps']}` تکرار و اجرای عمیق‌تر با `{ctx['bootstrap_deeper_reps']}` تکرار را ارائه می‌کند. مقادیر میانگین، اختلاف متوسط، و هم‌بستگی بین اجراها در آن آمده است تا پایداری نتایج نسبت به عمق بوت‌استرپ مستند شود.",
            "en": f"This table provides the summary comparison between the baseline bootstrap with `{ctx['bootstrap_reps']}` replicates and the deeper rerun with `{ctx['bootstrap_deeper_reps']}` replicates. Means, average differences, and cross-run correlations are reported to document sensitivity to bootstrap depth.",
        },
        "bootstrap_depth_sensitivity_station_comparison.csv": {
            "fa": f"این جدول مقایسه‌ی ایستگاه‌به‌ایستگاه بین دو عمق بوت‌استرپ را در یک قاب مشترک ارائه می‌کند و خروجی‌های `{ctx['bootstrap_reps']}` و `{ctx['bootstrap_deeper_reps']}` تکرار را کنار هم قرار می‌دهد. این فایل برای رهگیری اختلاف‌های موضعی در برآوردهای عدم‌قطعیت مفید است.",
            "en": f"This table presents the station-by-station comparison between the two bootstrap depths in a common frame, placing the `{ctx['bootstrap_reps']}`- and `{ctx['bootstrap_deeper_reps']}`-replicate outputs side by side. It is useful for tracing local differences in uncertainty estimates.",
        },
        "bootstrap_method_sensitivity_station_level.csv": {
            "fa": "این جدول حساسیت روش بوت‌استرپ را در سطح ایستگاه مستند می‌کند و نتایج اجرای پایه را با روش‌های جایگزین بوت‌استرپ در یک قاب مقایسه‌ای کنار هم می‌گذارد. هدف آن ثبت تفاوت‌های موضعی برآوردهای مبتنی بر بوت‌استرپ در ایستگاه‌ها و شاخص‌های مختلف است.",
            "en": "This table documents bootstrap-method sensitivity at the station level by placing the baseline results and the alternative bootstrap-method outputs into a common comparison frame. Its purpose is to preserve local differences in bootstrap-based estimates across stations and indices.",
        },
        "bootstrap_method_sensitivity_summary.csv": {
            "fa": "این جدول جمع‌بندی حساسیت به روش بوت‌استرپ را ارائه می‌کند و برای هر شاخص و هر معیار، اندازه‌ی اختلاف متوسط بین روش پایه و روش جایگزین را گزارش می‌نماید. این فایل برای نمایش اینکه انتخاب روش بوت‌استرپ تا چه حد بر خلاصه‌های کلیدی اثر می‌گذارد استفاده می‌شود.",
            "en": "This table provides the summary of bootstrap-method sensitivity, reporting the average magnitude of differences between the baseline method and the alternative bootstrap method for each index and metric. It is used to show how strongly the choice of bootstrap method affects the key summaries.",
        },
        "driver_analysis_summary.csv": {
            "fa": "این جدول نتایج تحلیل متغیرهای توضیحی را ذخیره می‌کند و برای هر شاخص، هر معیار روند، و هر پیش‌بین مکانی، ضرایب استانداردشده، `p-value`ها، `R^2` مدل، و هم‌بستگی اسپیرمن را گزارش می‌کند. این فایل برای مستندسازی این به‌کار می‌رود که آیا الگوهای روند با متغیرهایی مانند عرض جغرافیایی، طول جغرافیایی و ارتفاع ارتباط دارند یا نه.",
            "en": "This table stores the driver-analysis results and reports standardized coefficients, `p-values`, model `R^2`, and Spearman correlations for each index, trend metric, and spatial predictor. It documents whether the trend patterns are associated with variables such as latitude, longitude, and elevation.",
        },
        "interpolation_method_sensitivity_summary.csv": {
            "fa": "این جدول حساسیت به روش درون‌یابی را خلاصه می‌کند و برای هر شاخص، هر صدک، و هر جفت از روش‌های درون‌یابی، هم‌بستگی سطحی و `RMSE` بین سطوح تولیدشده را ثبت می‌نماید. این خروجی نشان می‌دهد که سطوح نمایشی نقشه‌ها تا چه حد به انتخاب روش درون‌یابی وابسته‌اند.",
            "en": "This table summarizes interpolation-method sensitivity by recording surface correlation and `RMSE` between pairs of interpolated surfaces for each index and quantile. It shows the extent to which the visual map surfaces depend on the interpolation choice.",
        },
        "regional_cluster_composites.csv": {
            "fa": "این جدول خلاصه‌های مرکب خوشه‌ای را نگه می‌دارد و برای هر شاخص، خوشه و معیار، آماره‌هایی مانند میانگین، میانه، انحراف معیار، کمینه و بیشینه را ارائه می‌کند. این فایل توصیف عددی هر ناحیه‌ی خوشه‌ای را به‌عنوان یک واحد منطقه‌ای فراهم می‌کند.",
            "en": "This table stores the cluster-composite summaries and reports statistics such as mean, median, standard deviation, minimum, and maximum for each index, cluster, and metric. It provides the numerical description of each cluster-defined regional unit.",
        },
        "regional_cluster_spatial_validation.csv": {
            "fa": "این جدول اعتبارسنجی فضایی خوشه‌ها را ثبت می‌کند و فاصله‌ی متوسط درون‌خوشه‌ای مشاهده‌شده را با توزیع حاصل از برچسب‌گذاری تصادفی مقایسه می‌نماید. بدین ترتیب، فایل حاضر نشان می‌دهد که آیا خوشه‌ها از نظر فضایی نیز فشرده‌تر از حالت تصادفی هستند یا نه.",
            "en": "This table records the spatial validation of the clusters by comparing the observed mean within-cluster distance against the permutation-based reference distribution. It therefore shows whether the clusters are spatially more compact than expected under random label assignment.",
        },
        "representative_station_selection.csv": {
            "fa": "این جدول ایستگاه‌های نماینده‌ی هر خوشه را ثبت می‌کند؛ یعنی ایستگاه‌هایی که در فضای ویژگی غربال‌شده نزدیک‌ترین اعضای مشاهده‌شده به مراکز خوشه بوده‌اند. این فایل مبنای شکل‌های مقایسه‌ی ایستگاه‌های نماینده در بخش تکمیلی است.",
            "en": "This table records the representative stations of each cluster, defined as the observed stations lying closest to the corresponding cluster centroids in the screened feature space. It forms the basis of the representative-station comparison figures in the supplementary material.",
        },
        "spatial_autocorrelation_moran.csv": {
            "fa": "این جدول نتایج خودهمبستگی فضایی را در قالب `Moran's I` برای میدان‌های شیب ایستگاهی ثبت می‌کند. برای هر شاخص و هر صدک، مقدار `Moran's I`، تعداد ایستگاه‌ها، و `p-value` مبتنی بر جایگشت ذخیره شده است تا ساختار فضایی روندها به‌طور رسمی مستند شود.",
            "en": "This table records the spatial-autocorrelation results in terms of `Moran's I` for the station-slope fields. For each index and quantile, it stores the `Moran's I` value, the number of stations, and the permutation-based `p-value`, thereby formally documenting the spatial structure of the trends.",
        },
        "station_significance_fdr.csv": {
            "fa": "این جدول نتایج معنی‌داری ایستگاهی را پس از کنترل چندآزمونی نگه می‌دارد. برای هر ایستگاه، شاخص و صدک، شیب، بازه‌های اطمینان، `p-value` تحلیلی، مقدار `q` و وضعیت ردشدن در روش `FDR` گزارش شده است.",
            "en": "This table stores the station-level significance results after multiple-testing control. For each station, index, and quantile, it reports the slope, confidence intervals, analytic `p-value`, `q` value, and the rejection status under the `FDR` procedure.",
        },
    }
    return descriptions.get(name, fallback_table_desc(name))


def fallback_figure_desc(name: str) -> dict[str, str]:
    return {
        "kind_fa": "شکل",
        "kind_en": "Figure",
        "fa": f"این شکل یکی از خروجی‌های گرافیکی مستقیم پایپ‌لاین است و فایل `{name}` را به‌عنوان بخشی از مستندات تکمیلی تحلیل ثبت می‌کند.",
        "en": f"This figure is a direct graphical output of the pipeline and documents `{name}` as part of the supplementary analytical record.",
    }


def describe_root_figure(figure_path: Path, ctx: dict) -> dict[str, str]:
    name = figure_path.name
    p = ctx["primary_delta"]
    if name == "data_coverage_by_station.png":
        return {
            "kind_fa": "نمودار پوشش داده",
            "kind_en": "Data-coverage bar chart",
            "fa": f"این شکل تعداد مقادیر سالانه‌ی معتبر موجود برای هر ایستگاه را در بازه‌ی `{ctx['analysis_years']}` نشان می‌دهد. مقادیر پس از ساخت شاخص‌های سالانه شمارش شده‌اند و هدف شکل، مستندسازی طول رکورد قابل‌استفاده و پوشش داده‌ی ایستگاه‌ها است.",
            "en": f"This figure shows the number of valid annual values available for each station during `{ctx['analysis_years']}`. The counts are derived after annual-index construction, and the purpose of the figure is to document record length and usable data coverage across stations.",
        }
    m = re.match(r"region_quantile_slopes_(.+)\.png", name)
    if m:
        idx_name = m.group(1)
        return {
            "kind_fa": "پروفایل ضرایب صدکی منطقه‌ای",
            "kind_en": "Regional quantile-coefficient profile",
            "fa": f"این شکل پروفایل شیب رگرسیون صدکی را برای میانگین سالانه‌ی همه‌ی ایستگاه‌ها در شاخص «{index_title_fa(idx_name)}» نمایش می‌دهد. ابتدا میانگین شبکه‌ای سالانه ساخته شده، سپس رگرسیون صدکی بر روی شبکه‌ی کامل صدک‌ها اجرا شده و شیب `OLS` نیز به‌عنوان مرجع میانگین در همان قاب آمده است.",
            "en": f"This figure displays the quantile-regression slope profile for the annual network-mean series of {index_title_en(idx_name)}. The annual station mean was first aggregated across the network, after which quantile regression was fitted across the full quantile grid, with the `OLS` slope shown as the mean-reference benchmark.",
        }
    if name == "ijoc_regional_quantile_panels.png":
        return {
            "kind_fa": "شکل چندپنلی پروفایل‌های صدکی منطقه‌ای",
            "kind_en": "Multi-panel regional quantile-profile figure",
            "fa": "این شکل چهار شاخص دمایی را در قالب یک مجموعه‌ی چندپنلی از پروفایل‌های صدکی منطقه‌ای کنار هم قرار می‌دهد. هر پنل مبتنی بر رگرسیون صدکی روی سری میانگین سالانه‌ی شبکه است و برای ارائه‌ی یک نمای توصیفی و مقاله‌ای از رفتار توزیعی شاخص‌ها طراحی شده است.",
            "en": "This figure places the four thermal indices into a multi-panel set of regional quantile profiles. Each panel is based on quantile regression fitted to the annual network-mean series and is intended to provide a descriptive, manuscript-style view of distributional trend structure.",
        }
    if name == "ijoc_study_area.png":
        return {
            "kind_fa": "نقشه‌ی منطقه‌ی مطالعه",
            "kind_en": "Study-area map",
            "fa": "این نقشه محدوده‌ی مطالعه و شبکه‌ی ایستگاه‌های هواشناسی را بر روی مرز ایران نمایش می‌دهد. رنگ‌آمیزی نقاط بر اساس ارتفاع انجام شده است تا توزیع مکانی ایستگاه‌ها و زمینه‌ی توپوگرافی مطالعه به‌صورت مقدماتی و رسمی معرفی شود.",
            "en": "This map presents the study domain and the meteorological station network over the national boundary of Iran. Station points are colored by elevation in order to introduce both the spatial station distribution and the topographic context of the study in a formal overview.",
        }
    m = re.match(r"station_focus_heatmap_(.+)\.png", name)
    if m:
        idx_name = m.group(1)
        return {
            "kind_fa": "هیت‌مپ امضای روند ایستگاهی",
            "kind_en": "Station trend-signature heatmap",
            "fa": f"این هیت‌مپ امضای روند ایستگاه‌ها را برای شاخص «{index_title_fa(idx_name)}» خلاصه می‌کند. ستون‌ها شامل شیب‌های صدکی کانونی `{ctx['focus_quantiles_text']}` و دلتای اصلی `{p}` هستند و سطرها ایستگاه‌ها را در قالبی مناسب برای مقایسه‌ی الگوهای نسبی بین ایستگاه‌ها نشان می‌دهند.",
            "en": f"This heatmap summarizes station trend signatures for {index_title_en(idx_name)}. The columns contain the focal quantile slopes `{ctx['focus_quantiles_text']}` and the main delta metric `{p}`, while the rows display stations in a format suitable for comparing relative pattern structure across the network.",
        }
    m = re.match(r"delta_uncertainty_(.+)\.png", name)
    if m:
        idx_name = m.group(1)
        return {
            "kind_fa": "نمودار عدم‌قطعیت دلتا",
            "kind_en": "Delta-uncertainty figure",
            "fa": f"این شکل میانگین بوت‌استرپی `{p}` و بازه‌ی اطمینان `95%` آن را برای شاخص «{index_title_fa(idx_name)}» در سطح ایستگاه نمایش می‌دهد. شکل حاضر صرفا ساختار برآورد مرکزی و دامنه‌ی عدم‌قطعیت را معرفی می‌کند و بر بوت‌استرپ اصلی با روش `{ctx['bootstrap_method']}` و `{ctx['bootstrap_reps']}` تکرار تکیه دارد.",
            "en": f"This figure displays the bootstrap mean of `{p}` together with its `95%` confidence interval for {index_title_en(idx_name)} at the station level. It is intended to describe the structure of the central estimate and uncertainty range, based on the main bootstrap using the `{ctx['bootstrap_method']}` method with `{ctx['bootstrap_reps']}` replicates.",
        }
    m = re.match(r"dendrogram_(.+)\.png", name)
    if m:
        idx_name = m.group(1)
        return {
            "kind_fa": "دندروگرام خوشه‌بندی پایه",
            "kind_en": "Baseline clustering dendrogram",
            "fa": f"این دندروگرام ساختار شباهت بین ایستگاه‌ها را برای شاخص «{index_title_fa(idx_name)}» در خوشه‌بندی پایه نشان می‌دهد. شاخه‌ها حاصل اجرای خوشه‌بندی `{ctx['cluster_algorithm']}` با `linkage={ctx['cluster_linkage']}` و `metric={ctx['cluster_metric']}` هستند و طول آن‌ها فاصله‌ی نسبی بین امضاهای روند را نمایش می‌دهد.",
            "en": f"This dendrogram shows the station-similarity structure for {index_title_en(idx_name)} under the baseline clustering configuration. The branches arise from `{ctx['cluster_algorithm']}` clustering with `linkage={ctx['cluster_linkage']}` and `metric={ctx['cluster_metric']}`, and branch length represents relative separation among station trend signatures.",
        }
    if name == "ijoc_data_quality_homogeneity.png":
        return {
            "kind_fa": "شکل ترکیبی کیفیت داده و همگنی",
            "kind_en": "Combined data-quality and homogeneity figure",
            "fa": "این شکل در یک قاب دوپنلی، کامل‌بودن داده‌های دمایی و خلاصه‌ی نتایج آزمون‌های همگنی را نشان می‌دهد. پنل نخست توزیع کامل‌بودن ایستگاه‌ها را نمایش می‌دهد و پنل دوم شمار ایستگاه‌های پرچم‌خورده در آزمون‌های همگنی روندزدایی‌شده را خلاصه می‌کند.",
            "en": "This figure combines data completeness and homogeneity diagnostics within a two-panel layout. The first panel shows the distribution of station-level completeness, whereas the second summarizes the number of stations flagged by the detrended homogeneity tests.",
        }
    if name == "ijoc_homogeneity_sensitivity.png":
        return {
            "kind_fa": "هیت‌مپ حساسیت به حذف ایستگاه‌های پرچم‌خورده",
            "kind_en": "Sensitivity heatmap for flagged-station exclusion",
            "fa": f"این شکل اختلاف خلاصه‌های منطقه‌ای را بین اجرای کامل و اجرای بدون ایستگاه‌های پرچم‌خورده در آزمون همگنی نشان می‌دهد. مقادیر برای شیب‌های `OLS`، شیب‌های صدکی کانونی و `{p}` گزارش شده‌اند تا اثر این انتخاب روش‌شناختی در قالبی مکمل و فشرده ثبت شود.",
            "en": f"This figure shows the differences in regional summaries between the full analysis and the rerun excluding homogeneity-flagged stations. Values are reported for `OLS` slopes, focal quantile slopes, and `{p}` in order to document the effect of this methodological choice in a concise supplementary form.",
        }
    if name == "ijoc_bootstrap_depth_sensitivity.png":
        return {
            "kind_fa": "شکل حساسیت عمق بوت‌استرپ",
            "kind_en": "Bootstrap-depth sensitivity figure",
            "fa": f"این شکل خلاصه‌های منطقه‌ای بوت‌استرپ را بین اجرای پایه با `{ctx['bootstrap_reps']}` تکرار و اجرای عمیق‌تر با `{ctx['bootstrap_deeper_reps']}` تکرار مقایسه می‌کند. پنل‌ها برای مستندسازی پایداری نتایج نسبت به افزایش تعداد تکرارهای بوت‌استرپ طراحی شده‌اند.",
            "en": f"This figure compares regional bootstrap summaries between the baseline run with `{ctx['bootstrap_reps']}` replicates and the deeper rerun with `{ctx['bootstrap_deeper_reps']}` replicates. The panels are designed to document the stability of the results with respect to increased bootstrap depth.",
        }
    if name == "ijoc_alternative_clustering_sensitivity.png":
        return {
            "kind_fa": "هیت‌مپ حساسیت روش‌های خوشه‌بندی",
            "kind_en": "Alternative-clustering sensitivity heatmap",
            "fa": "این شکل میزان توافق خوشه‌بندی پایه با چند روش جایگزین را در قالب هیت‌مپ `Adjusted Rand Index` نمایش می‌دهد. هر سلول رابطه‌ی بین یک شاخص و یک روش خوشه‌بندی جایگزین را خلاصه می‌کند و برای مستندسازی پایداری منطقه‌بندی استفاده می‌شود.",
            "en": "This figure presents the agreement between the baseline clustering solution and several alternative methods as an `Adjusted Rand Index` heatmap. Each cell summarizes the relationship between one index and one alternative clustering specification, thereby documenting regionalization stability.",
        }
    if name == "ijoc_main_delta1_maps.png":
        return {
            "kind_fa": "نقشه‌ی چندپنلی دلتای اصلی",
            "kind_en": "Multi-panel main-delta map",
            "fa": f"این شکل توزیع فضایی `{p}` را برای هر چهار شاخص در قالب یک مجموعه‌ی چندپنلی از نقشه‌های ایستگاهی نشان می‌دهد. رنگ هر نقطه مستقیما از مقدار ثبت‌شده در جدول ویژگی‌ها گرفته شده و این خروجی به‌عنوان نمای اصلی فضایی معیار دلتا در اسناد تکمیلی عمل می‌کند.",
            "en": f"This figure shows the spatial distribution of `{p}` for all four indices in a multi-panel set of station maps. Point colors are taken directly from the values stored in the feature table, making this output the principal spatial presentation of the delta metric in the supplementary record.",
        }
    if name == "ijoc_split_period_comparison.png":
        return {
            "kind_fa": "شکل مقایسه‌ی دو زیر‌دوره",
            "kind_en": "Split-period comparison figure",
            "fa": f"این شکل شیب‌های صدکی کانونی و `OLS` را بین دو زیر‌دوره‌ی `{ctx['split_1']}` و `{ctx['split_2']}` مقایسه می‌کند. تقسیم بازه‌ی تحلیل به دو نیمه برای مستندسازی دیداری تفاوت ساختار روند بین بخش‌های ابتدایی و انتهایی دوره به‌کار رفته است.",
            "en": f"This figure compares focal quantile slopes and `OLS` slopes between the two sub-periods `{ctx['split_1']}` and `{ctx['split_2']}`. The analysis period is split into halves to provide a visual supplementary record of how trend structure differs between the earlier and later parts of the record.",
        }
    if name == "ijoc_robustness_synthesis.png":
        return {
            "kind_fa": "شکل ترکیبی پایداری و حساسیت",
            "kind_en": "Robustness-synthesis figure",
            "fa": "این شکل چندپنلی، خروجی‌های اصلی حساسیت روش‌شناختی و پایداری خوشه‌بندیِ موجود در اجرای فعلی را در یک قاب واحد کنار هم می‌گذارد. بسته به این‌که برای هر مؤلفه داده‌ی مقایسه‌ای تولید شده باشد، پنل‌های مربوط به دوره‌ی مرجع، روش بوت‌استرپ، روش درون‌یابی و پایداری خوشه‌بندی به‌صورت پویا نمایش داده می‌شوند.",
            "en": "This multi-panel figure assembles the available methodological-sensitivity and clustering-robustness outputs from the current run into a single layout. Depending on which comparison products were generated, panels for reference period, bootstrap method, interpolation method, and clustering robustness are shown dynamically.",
        }
    m = re.match(r"map_(.+)_cluster\.png", name)
    if m:
        idx_name = m.group(1)
        return {
            "kind_fa": "نقشه‌ی خوشه‌های ایستگاهی",
            "kind_en": "Station-cluster map",
            "fa": f"این نقشه برچسب خوشه‌ی هر ایستگاه را برای شاخص «{index_title_fa(idx_name)}» نمایش می‌دهد. رنگ‌ها از اجرای پایه‌ی خوشه‌بندی گرفته شده‌اند و نقش شکل، نمایش پراکنش فضایی گروه‌های رفتاری حاصل از منطقه‌بندی است.",
            "en": f"This map displays the cluster label of each station for {index_title_en(idx_name)}. The colors are taken from the baseline clustering solution, and the figure is intended to show the spatial distribution of the regionalized behavioral groups.",
        }
    m = re.match(r"map_(.+)_Delta1\.png", name)
    if m:
        idx_name = m.group(1)
        return {
            "kind_fa": "نقشه‌ی ایستگاهی دلتا",
            "kind_en": "Station-based delta map",
            "fa": f"این نقشه مقدار ایستگاهی `{p}` را برای شاخص «{index_title_fa(idx_name)}» به‌صورت پراکنش نقطه‌ای نمایش می‌دهد. هیچ سطح درون‌یابی‌شده‌ای در این فایل به‌کار نرفته و هر نقطه مستقیما بیانگر مقدار ایستگاه متناظر است.",
            "en": f"This map shows the station-level values of `{p}` for {index_title_en(idx_name)} as a point-based spatial display. No interpolated surface is used in this file; each point directly represents the corresponding station value.",
        }
    m = re.match(r"map_(.+)_slope_(0\.\d+)\.png", name)
    if m:
        idx_name = m.group(1)
        tau = tau_label(m.group(2))
        return {
            "kind_fa": "نقشه‌ی ایستگاهی شیب صدکی",
            "kind_en": "Station-based quantile-slope map",
            "fa": f"این نقشه شیب رگرسیون صدکی در `τ={tau}` را برای شاخص «{index_title_fa(idx_name)}» در سطح ایستگاه‌ها نشان می‌دهد. رنگ هر نقطه از برآورد صدکی همان ایستگاه گرفته شده است و بنابراین شکل، نمایش مستقیم الگوی فضایی مقادیر مشاهده‌ای ایستگاهی محسوب می‌شود.",
            "en": f"This map presents the quantile-regression slope at `τ={tau}` for {index_title_en(idx_name)} at the station level. Each point is colored by the station-specific quantile estimate, so the figure constitutes a direct spatial display of observed station-based values.",
        }
    return fallback_figure_desc(name)


def describe_nested_figure(figure_path: Path, ctx: dict) -> dict[str, str] | None:
    parts = figure_path.relative_to(FIGURES).parts
    name = figure_path.name
    if parts[0] == "ijoc_station_comparisons":
        if name == "main_figure_representative_stations.png":
            return {
                "kind_fa": "شکل چندپنلی ایستگاه‌های نماینده",
                "kind_en": "Multi-panel representative-station figure",
                "fa": "این شکل ایستگاه‌های نماینده‌ی خوشه‌ها را برای هر چهار شاخص در یک قاب چندپنلی گردآوری می‌کند. ایستگاه‌های نماینده به‌صورت الگوریتمی و بر پایه‌ی نزدیک‌ترین عضو مشاهده‌شده به مرکز خوشه در فضای ویژگی انتخاب شده‌اند.",
                "en": "This figure assembles the representative stations of the clusters for all four indices within a multi-panel layout. Representative stations were selected algorithmically as the observed members closest to the corresponding cluster centroids in feature space.",
            }
        m = re.match(r"station_comparison_(.+)\.png", name)
        if m:
            idx_name = m.group(1)
            return {
                "kind_fa": "نمودار مقایسه‌ی ایستگاه‌های نماینده",
                "kind_en": "Representative-station comparison plot",
                "fa": f"این شکل پروفایل‌های کامل شیب رگرسیون صدکی را برای چند ایستگاه نماینده‌ی خوشه‌ها در شاخص «{index_title_fa(idx_name)}» کنار هم قرار می‌دهد. هدف آن معرفی شکل‌های شاخص رفتار خوشه‌ای در قالب نمونه‌های ایستگاهی قابل‌ردیابی است.",
                "en": f"This figure places the full quantile-regression slope profiles of several cluster-representative stations side by side for {index_title_en(idx_name)}. Its purpose is to present recognizable station-based examples of the cluster-defined behavior patterns.",
            }
    if parts[0] == "paper2_figure3_quantile_maps":
        m = re.match(r"figure3_tau_(0\.\d+)\.png", name)
        if m:
            tau = tau_label(m.group(1))
            return {
                "kind_fa": "نقشه‌ی درون‌یابی‌شده‌ی صدکی",
                "kind_en": "Interpolated quantile map",
                "fa": f"این شکل نقشه‌های شیب در `τ={tau}` را برای هر چهار شاخص نشان می‌دهد. پس از برآورد شیب صدکی در هر ایستگاه، یک سطح درون‌یابی‌شده با روش `{ctx['interp_method']}` و پارامتر هموارسازی `{ctx['interp_smooth']}` صرفا برای نمایش بصری ترسیم شده و نقاط ایستگاهی به‌عنوان شواهد اصلی روی آن قرار گرفته‌اند.",
                "en": f"This figure displays the slope maps at `τ={tau}` for all four indices. After estimating the station-level quantile slope, an interpolated surface was drawn using the `{ctx['interp_method']}` method with smoothing parameter `{ctx['interp_smooth']}` for visual orientation only, while the station points remain the primary evidence.",
            }
    if parts[0] == "paper2_station_figures":
        if ctx["paper2_fig1_sample"] and figure_path == ctx["paper2_fig1_sample"]:
            return {
                "kind_fa": "نمونه‌ی شکل سری زمانی ایستگاهی",
                "kind_en": "Sample station time-series figure",
                "fa": f"این فایل یک نمونه از مجموعه‌ی `station_figures/figure1_timeseries` است و برای یک ایستگاه، سری سالانه‌ی چهار شاخص را به‌همراه خطوط رگرسیون صدکی در صدک‌های `{ctx['focus_quantiles_text']}` و خط `OLS` نمایش می‌دهد. همین ساختار برای سایر ایستگاه‌ها نیز با قالب یکسان تکرار شده است.",
                "en": f"This file is a sample from the `station_figures/figure1_timeseries` set and shows the annual series of the four indices for one station together with quantile-regression lines at `{ctx['focus_quantiles_text']}` and the `OLS` line. The same structure is repeated for the remaining stations in an identical format.",
            }
        if ctx["paper2_fig2_sample"] and figure_path == ctx["paper2_fig2_sample"]:
            return {
                "kind_fa": "نمونه‌ی شکل پروفایل ضرایب صدکی ایستگاهی",
                "kind_en": "Sample station quantile-profile figure",
                "fa": "این فایل یک نمونه از مجموعه‌ی `station_figures/figure2_quantile_coefficients` است و برای یک ایستگاه، پروفایل کامل شیب رگرسیون صدکی را در سراسر شبکه‌ی صدک‌ها برای هر چهار شاخص نشان می‌دهد. بقیه‌ی فایل‌های این مجموعه نیز با همین قالب و تنها با داده‌ی ایستگاه‌های دیگر تولید شده‌اند.",
                "en": "This file is a sample from the `station_figures/figure2_quantile_coefficients` set and shows the full quantile-regression slope profile across the quantile grid for one station and all four indices. The remaining files in the set follow the same format and differ only in station-specific data.",
            }
        return None
    if parts[0] == "paper1_station_figures":
        if ctx["paper1_fig4_sample"] and figure_path == ctx["paper1_fig4_sample"]:
            return {
                "kind_fa": "نمونه‌ی شکل توزیع بوت‌استرپ ایستگاهی",
                "kind_en": "Sample station bootstrap-distribution figure",
                "fa": "این فایل یک نمونه از مجموعه‌ی `figure4_bootstrap_distributions` است و توزیع بوت‌استرپی شیب‌ها را برای یک ایستگاه، چند صدک کانونی، و چهار شاخص نمایش می‌دهد. در این قالب، هیستوگرام، چگالی تقریبی، برآورد نقطه‌ای، و بازه‌ی اطمینان بوت‌استرپی در یک قاب واحد ارائه شده‌اند.",
                "en": "This file is a sample from the `figure4_bootstrap_distributions` set and displays the bootstrap distributions of slopes for one station, several focal quantiles, and all four indices. The layout combines histogram-based density, approximate smoothing, the point estimate, and the bootstrap confidence interval within a single frame.",
            }
        return None
    if parts[0] == "paper1_quantile_dendrograms":
        m = re.match(r"figure_(\d+)_tau_(0\.\d+)\.png", name)
        if m:
            tau = tau_label(m.group(2))
            return {
                "kind_fa": "دندروگرام تک‌صدکی",
                "kind_en": "Single-quantile dendrogram",
                "fa": f"این شکل دندروگرام‌های میانگین‌پیوندی را برای شیب‌های ایستگاهی در `τ={tau}` ارائه می‌کند. در هر پنل، خوشه‌بندی صرفا بر اساس فاصله‌ی بین مقادیر شیب همان صدک انجام شده است تا ساختار شباهت تک‌صدکی به‌صورت مستقل نمایش داده شود.",
                "en": f"This figure presents average-linkage dendrograms for station slopes at `τ={tau}`. Within each panel, clustering is performed solely on the distances among slope values at that quantile so that the single-quantile similarity structure can be displayed independently.",
            }
    if parts[0] == "advanced_spatial_inference":
        if name == "raw_vs_fdr_counts.png":
            return {
                "kind_fa": "هیت‌مپ شمارش معنی‌داری خام و FDR",
                "kind_en": "Raw-vs-FDR significance heatmap",
                "fa": "این شکل تعداد نتایج معنی‌دار ایستگاهی را پیش و پس از کنترل چندآزمونی `FDR` مقایسه می‌کند. دو پنل آن نشان می‌دهند که در هر شاخص و صدک، چه تعداد نتیجه در غربال‌گری خام معنی‌دار بوده و چه تعداد پس از تصحیح `FDR` باقی مانده است.",
                "en": "This figure compares the number of significant station-level results before and after `FDR` multiple-testing control. Its two panels show, for each index and quantile, how many results are significant under raw screening and how many remain after `FDR` correction.",
            }
        if name == "moran_i_heatmap.png":
            return {
                "kind_fa": "هیت‌مپ Moran's I",
                "kind_en": "Moran's I heatmap",
                "fa": "این شکل مقادیر `Moran's I` را برای میدان‌های شیب ایستگاهی در شاخص‌ها و صدک‌های مختلف نمایش می‌دهد. نقش آن ارائه‌ی یک جمع‌بندی بصری از شدت و جهت خودهمبستگی فضایی در نتایج روند است.",
                "en": "This figure displays the `Moran's I` values for the station-slope fields across indices and quantiles. Its role is to provide a visual summary of the strength and direction of spatial autocorrelation in the trend results.",
            }
    if parts[0] == "advanced_method_sensitivity":
        if name == "bootstrap_method_sensitivity.png":
            return {
                "kind_fa": "شکل حساسیت روش بوت‌استرپ",
                "kind_en": "Bootstrap-method sensitivity figure",
                "fa": "این شکل اندازه‌ی اختلاف میان روش بوت‌استرپ پایه و روش جایگزین را برای معیارهای اصلی در شاخص‌های مختلف خلاصه می‌کند. این خروجی برای ارائه‌ی سریع میزان حساسیت خلاصه‌های عدم‌قطعیت نسبت به انتخاب روش بوت‌استرپ طراحی شده است.",
                "en": "This figure summarizes the magnitude of the differences between the baseline bootstrap method and the alternative method for the main metrics across the indices. It is designed to provide a rapid visual account of uncertainty sensitivity to bootstrap-method choice.",
            }
        m = re.match(r"interpolation_sensitivity_tau_(0\.\d+)\.png", name)
        if m:
            tau = tau_label(m.group(1))
            return {
                "kind_fa": "شکل حساسیت روش درون‌یابی",
                "kind_en": "Interpolation-sensitivity figure",
                "fa": f"این شکل سطوح درون‌یابی‌شده‌ی تولیدشده با چند روش را در `τ={tau}` کنار هم قرار می‌دهد تا تفاوت‌های دیداری بین روش‌های درون‌یابی مستند شود. این فایل تنها برای مقایسه‌ی روش‌های نمایش فضایی ساخته شده و جایگزین شواهد ایستگاهی نیست.",
                "en": f"This figure places the interpolated surfaces generated by several methods side by side at `τ={tau}` in order to document the visual differences among interpolation schemes. It is intended solely for comparing spatial-display methods and not as a substitute for the station-based evidence.",
            }
    if parts[0] == "advanced_driver_analysis":
        if name == "driver_effects_delta1_heatmap.png":
            return {
                "kind_fa": "هیت‌مپ اثر متغیرهای توضیحی بر Delta1",
                "kind_en": "Driver-effects heatmap for Delta1",
                "fa": "این شکل ضرایب استانداردشده‌ی متغیرهای توضیحی را برای معیار `Delta1` در شاخص‌های مختلف به‌صورت هیت‌مپ نشان می‌دهد. هدف آن ارائه‌ی یک نمای فشرده از جهت و بزرگی نسبی اثر پیش‌بین‌های مکانی بر معیار دلتا است.",
                "en": "This figure shows the standardized driver coefficients for the `Delta1` metric across the different indices in heatmap form. Its purpose is to provide a compact view of the direction and relative magnitude of the spatial predictors’ effects on the delta metric.",
            }
        m = re.match(r"driver_scatter_(.+)\.png", name)
        if m:
            idx_name = m.group(1)
            return {
                "kind_fa": "نمودار پراکنش متغیرهای توضیحی",
                "kind_en": "Driver scatterplot figure",
                "fa": f"این شکل روابط پراکنشی بین پیش‌بین‌های مکانی و معیار `Delta1` را برای شاخص «{index_title_fa(idx_name)}» نشان می‌دهد. هر پنل یکی از متغیرهای توضیحی را نمایش می‌دهد و خط برازش صرفا برای نمایش جهت کلی رابطه افزوده شده است.",
                "en": f"This figure shows scatterplot relationships between the spatial predictors and the `Delta1` metric for {index_title_en(idx_name)}. Each panel represents one predictor, and the fitted line is included only to visualize the overall direction of association.",
            }
    if parts[0] == "advanced_regionalization":
        m = re.match(r"cluster_composite_heatmap_(.+)\.png", name)
        if m:
            idx_name = m.group(1)
            return {
                "kind_fa": "هیت‌مپ مرکب خوشه‌ای",
                "kind_en": "Cluster-composite heatmap",
                "fa": f"این شکل میانه‌ی معیارهای منتخب را برای خوشه‌های مختلف در شاخص «{index_title_fa(idx_name)}» نشان می‌دهد. هیت‌مپ حاضر یک نمای جمع‌بندی از تفاوت ساختار کمی خوشه‌ها در فضای معیارهای روند ارائه می‌کند.",
                "en": f"This figure shows the median values of the selected metrics across clusters for {index_title_en(idx_name)}. The heatmap provides a summary view of the quantitative differences among clusters in trend-metric space.",
            }
        m = re.match(r"cluster_delta1_boxplot_(.+)\.png", name)
        if m:
            idx_name = m.group(1)
            return {
                "kind_fa": "جعبه‌نمودار Delta1 در خوشه‌ها",
                "kind_en": "Cluster-wise Delta1 boxplot",
                "fa": f"این شکل توزیع `Delta1` را به تفکیک خوشه برای شاخص «{index_title_fa(idx_name)}» نمایش می‌دهد. هر جعبه‌نمودار دامنه‌ی تغییرپذیری درون‌خوشه‌ای را نشان می‌دهد و برای مقایسه‌ی دیداری تفاوت خوشه‌ها در معیار دلتا به‌کار می‌رود.",
                "en": f"This figure displays the distribution of `Delta1` by cluster for {index_title_en(idx_name)}. Each boxplot shows the within-cluster variability and is used for visually comparing cluster differences in the delta metric.",
            }
        if name == "cluster_spatial_validation.png":
            return {
                "kind_fa": "شکل اعتبارسنجی فضایی خوشه‌ها",
                "kind_en": "Cluster spatial-validation figure",
                "fa": "این شکل فاصله‌ی متوسط درون‌خوشه‌ای مشاهده‌شده را با میانگین فاصله‌ی حاصل از جایگشت تصادفی برای هر شاخص مقایسه می‌کند. یادداشت‌های بالای ستون‌ها نیز `p-value` آزمون فشردگی فضایی را نشان می‌دهند.",
                "en": "This figure compares the observed mean within-cluster distance with the permutation-based mean distance for each index. The annotations above the bars also report the `p-value` of the spatial-compactness test.",
            }
    return None


def generate() -> None:
    cfg = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    ctx = build_context(cfg)

    if DOCS.exists():
        for old in sorted(DOCS.rglob("*.md"), reverse=True):
            old.unlink()

    figure_docs: list[tuple[str, Path]] = []
    table_docs: list[tuple[str, Path]] = []

    for table_path in sorted(TABLES.glob("*.csv")):
        doc_path = DOCS / "tables" / f"{table_path.stem}.md"
        table_doc(doc_path, table_path, describe_table(table_path, ctx))
        table_docs.append((table_path.name, doc_path))

    for figure_path in sorted(FIGURES.glob("*.png")):
        doc_path = DOCS / "figures" / f"{figure_path.stem}.md"
        figure_doc(doc_path, figure_path, describe_root_figure(figure_path, ctx))
        figure_docs.append((figure_path.name, doc_path))

    for figure_path in sorted(FIGURES.rglob("*.png")):
        if figure_path.parent == FIGURES:
            continue
        desc = describe_nested_figure(figure_path, ctx)
        if not desc:
            continue
        doc_path = DOCS / "figures" / figure_path.relative_to(FIGURES).parent / f"{figure_path.stem}.md"
        figure_doc(doc_path, figure_path, desc)
        figure_docs.append((figure_path.relative_to(FIGURES).as_posix(), doc_path))

    figure_docs.sort(key=lambda x: x[0])
    table_docs.sort(key=lambda x: x[0])

    lines = [
        "# مستندات تکمیلی خروجی‌های پایپ‌لاین / Supplementary Output Documentation",
        "",
        "## فارسی",
        "",
        f"این پوشه مجموعه‌ای از یادداشت‌های توصیفی دو‌زبانه برای خروجی‌های تولیدشده در اجرای فعلی در بازه‌ی `{ctx['analysis_years']}` را فراهم می‌کند. هدف این مستندات، معرفی ماهیت هر فایل، جایگاه آن در گردش‌کار، و روش کلی تولید آن است؛ بنابراین متن‌ها عمدا از تحلیل نتایج پرهیز می‌کنند و برای استفاده در بخش Supplementary/Appendix تنظیم شده‌اند.",
        "",
        "## English",
        "",
        f"This folder provides a bilingual set of descriptive notes for the outputs generated in the current run over `{ctx['analysis_years']}`. The purpose of these documents is to identify the role of each file, its position within the workflow, and the general method used to produce it; accordingly, the text intentionally avoids interpretation and is framed for Supplementary/Appendix use.",
        "",
        "## Note / نکته",
        "",
        "- For `station_figures`, only one representative example from each figure family is documented because the remaining files follow the same layout and differ only in station-specific data. / برای `station_figures` فقط یک نمونه‌ی نماینده از هر خانواده‌ی شکل مستندسازی شده است، زیرا بقیه‌ی فایل‌ها همان قالب را با داده‌ی ایستگاه‌های دیگر تکرار می‌کنند.",
        "",
        "## Figures And Maps / شکل‌ها و نقشه‌ها",
        "",
    ]
    for label, doc_path in figure_docs:
        lines.append(f"- [{label}]({make_rel(DOCS / 'README.md', doc_path)})")

    lines += [
        "",
        "## Tables / جدول‌ها",
        "",
    ]
    for label, doc_path in table_docs:
        lines.append(f"- [{label}]({make_rel(DOCS / 'README.md', doc_path)})")

    lines += [
        "",
        "## Auxiliary Outputs / خروجی‌های کمکی",
        "",
        f"- [REPORT.md]({make_rel(DOCS / 'README.md', OUTPUTS / 'REPORT.md')})",
        f"- [run_metadata.json]({make_rel(DOCS / 'README.md', OUTPUTS / 'run_metadata.json')})",
        f"- [run_status.txt]({make_rel(DOCS / 'README.md', OUTPUTS / 'run_status.txt')})",
        f"- [run_status_summary.txt]({make_rel(DOCS / 'README.md', OUTPUTS / 'run_status_summary.txt')})",
        f"- [run_status_detail.txt]({make_rel(DOCS / 'README.md', OUTPUTS / 'run_status_detail.txt')})",
    ]
    write_doc(DOCS / "README.md", "\n".join(lines))


if __name__ == "__main__":
    generate()
