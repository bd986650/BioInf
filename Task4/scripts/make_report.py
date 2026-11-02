#!/usr/bin/env python3
"""Generate the same concise Russian report as Markdown and a two-page PDF.

Only this optional presentation step requires reportlab; data processing uses
the Python standard library. Font paths can be specified on any platform.
"""
import argparse
import csv
import json
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[1]


def contents():
    with (ROOT/'results/annotations.tsv').open(encoding='utf-8') as f:
        rows = {r['original_id']: r for r in csv.DictReader(f, delimiter='\t')}
    stats = json.loads((ROOT/'results/statistics.json').read_text())
    def af(ident):
        r = rows[ident]
        return (f"{float(r['gnomad_genome_af']):.6g} "
                f"(AC/AN={r['gnomad_genome_ac']}/{r['gnomad_genome_an']})")
    return [
      ('title', 'Клиническая аннотация генетических вариантов'),
      ('subtitle', 'Домашнее задание 4 | Данил Белов | Вариант 6 | 02.11.2026'),
      ('heading', 'Исходные данные и статистика'),
      ('paragraph', 'Проанализирован student_06.vcf: человек, сборка GRCh38, VCFv4.2. '
       'Все записи - биаллельные SNV. CHROM - хромосома, POS - позиция с нумерацией от 1, '
       'ID - идентификатор, REF/ALT - референсный/альтернативный аллель. QUAL - качество, '
       'FILTER - фильтры, INFO - дополнительные аннотации; здесь эти три поля равны «.», '
       'то есть не заданы. Генотипов и данных пациента нет.'),
      ('table', [('Показатель', 'Результат'), ('Число вариантов', '50'),
       ('Хромосомы 21 / 22', '11 / 39'), ('Транзиции / трансверсии', '41 / 9; Ti/Tv = 4,56'),
       ('Максимальное последствие VEP', '32 missense; 18 synonymous'),
       ('ClinVar: запись есть / не найдена', '19 / 31')]),
      ('heading', 'Методика и настройки'),
      ('paragraph', 'Аннотация получена 02.11.2026 через публичные API: Ensembl VEP release 116 '
       '(homo_sapiens, GRCh38), ClinVar NCBI E-utilities, gnomAD GraphQL (dataset=gnomad_r4, GRCh38). '
       'VEP: canonical, mane, hgvs, protein, numbers, transcript_version, variant_class = 1; '
       'все транскрипты сохранены, pick и per_gene не применялись [1].'),
      ('paragraph', 'В итоговую строку переносится most_severe_consequence VEP. Среди транскриптов '
       'с этим последствием приоритет имеют MANE Select, затем canonical, затем protein_coding; '
       'последняя развязка - transcript_id по алфавиту. Это явное правило свёртки, а не клиническая '
       'оценка. У 9 вариантов максимальное последствие отличается от MANE Select; альтернативные '
       'последствия сохранены в таблицах.'),
      ('paragraph', 'ClinVar искали по rsID и независимо по хромосоме/позиции GRCh38 [2]. '
       'Из 24 выгруженных записей 19 совпали также по REF и ALT; 5 других изменений исключены. '
       'В сводной таблице приведена агрегированная germline-классификация VCV: 8 Benign, '
       '4 Benign/Likely benign, 2 Likely benign, 3 VUS и 2 конфликта. '
       'Для gnomAD приведены AF, AC, AN и фильтры отдельно для геномов и экзомов. '
       'Основная AF - genome при PASS и AN>0, иначе exome при тех же условиях; когорты не складываются. '
       'NA означает отсутствие данных, а не нулевую частоту.'),
      ('heading', '1. rs114053718 - SYNJ1: неопределённая значимость'),
      ('paragraph', '21:32656885 A>G; MANE NM_203446.3 / ENST00000674351.1, '
       'c.2597T>C, p.Ile866Thr: missense_variant. Различие REF/ALT в геномной и транскриптной '
       'записи объясняется отрицательной цепью гена. SIFT=0 (deleterious_low_confidence), '
       'PolyPhen=0.477 (possibly_damaging) - только вычислительные предсказания. '
       f'gnomAD genome AF={af("rs114053718")}; вариант редкий. '
       'ClinVar VCV000544554.8: Uncertain significance, три подателя без конфликта; '
       'последняя оценка 24.06.2024 [3]. Редкость и замена аминокислоты не доказывают патогенность.'),
      ('pagebreak', ''),
      ('heading', '2. rs116555717 - PDGFB: конфликт интерпретаций'),
      ('paragraph', '22:39231644 T>C; MANE NM_002608.4 / ENST00000331163.11, '
       'c.434A>G, p.Gln145Arg: missense_variant. SIFT=0.60 (tolerated), '
       'PolyPhen=0.039 (benign). '
       f'gnomAD genome AF={af("rs116555717")}. '
       'ClinVar VCV003305447.2: Conflicting classifications of pathogenicity - '
       'Uncertain significance (1) и Likely benign (1); последняя оценка 14.12.2025 [4]. '
       'Следовательно, согласованной клинической оценки нет. Запись VCV002585196 относится '
       'к T>G / p.Gln145Pro в той же позиции и здесь не используется. '
       'Предсказания SIFT/PolyPhen не разрешают конфликт ClinVar.'),
      ('heading', '3. rs1135638 - MRPL39: зависимость от транскрипта'),
      ('paragraph', '21:25592836 G>A. Максимальное последствие - missense_variant '
       'на ENST00000985908.1: c.833C>T, p.Ala278Val; биотип nonsense_mediated_decay, '
       'дополнительное последствие NMD_transcript_variant. Это аннотация уже существующего '
       'NMD-транскрипта, а не доказательство того, что SNV вызывает NMD. '
       'На MANE Select NM_017446.4 / ENST00000352957.9 вариант синонимический: '
       'c.897C>T, p.Gly299=. '
       f'gnomAD genome AF={af("rs1135638")}, то есть около 78,3% альтернативных аллелей [5]. '
       'Точная запись ClinVar не найдена; клиническая интерпретация отсутствует. '
       'Высокая частота ограничивает гипотезу высокопенетрантного редкого заболевания, '
       'но не заменяет клиническую классификацию.'),
      ('heading', 'Заключение и ограничения'),
      ('paragraph', 'Наиболее показательны SYNJ1 с VUS, PDGFB с конфликтом оценок и MRPL39 '
       'с различием между максимальным последствием и MANE. В наборе нет совпавших вариантов '
       'с агрегированной классификацией Pathogenic/Likely pathogenic; это не доказывает '
       'отсутствие генетического риска. Missense описывает молекулярную замену, а ClinVar - '
       'предоставленную подателями клиническую оценку; одно не определяет другое [6].'),
      ('paragraph', 'Ограничения: учебная выборка мала и неслучайна; Ti/Tv здесь не является '
       'оценкой качества всего генома. Нет генотипа, зиготности, покрытия, фенотипа, семейных '
       'и функциональных данных. Частоты gnomAD зависят от когорты, происхождения участников, '
       'покрытия и фильтров; редкость сама по себе не означает патогенность. Отсутствие ClinVar '
       'не означает доброкачественность; конфликты, давность оценок и статус рецензирования '
       'сохранены. Выбор транскрипта меняет аннотацию; API и базы обновляются, поэтому исходные '
       'ответы, время запросов, настройки и SHA-256 сохранены. Анализ представляет учебную '
       'аннотацию без индивидуальной диагностики.'),
      ('heading', 'Источники и результаты'),
      ('small', '[1] Ensembl REST VEP: https://rest.ensembl.org/documentation/info/vep_region_post\n'
       '[2] Поиск ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/docs/help/\n'
       '[3] SYNJ1: https://www.ncbi.nlm.nih.gov/clinvar/variation/544554/\n'
       '[4] PDGFB: https://www.ncbi.nlm.nih.gov/clinvar/variation/3305447/\n'
       '[5] gnomAD: https://gnomad.broadinstitute.org/variant/21-25592836-G-A?dataset=gnomad_r4\n'
       '[6] Интерпретации ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/'),
      ('small', 'Полная таблица: results/annotations.tsv (50 строк); отдельные выгрузки - '
       'results/vep_transcripts.tsv, clinvar.tsv, gnomad.tsv; исходные ответы - raw/. '
       f'Сохранено {stats["transcript_row_count"]} транскриптных аннотаций. Настройки: results/annotation_settings.json.')
    ]


def find_font(explicit, candidates):
    for candidate in ([explicit] if explicit else candidates):
        if candidate and Path(candidate).exists():
            return candidate
    raise FileNotFoundError('Specify --font-regular/--font-bold with a Cyrillic TrueType font')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--font-regular')
    parser.add_argument('--font-bold')
    args = parser.parse_args()
    blocks = contents()
    markdown = []
    for kind, body in blocks:
        if kind == 'pagebreak':
            continue
        if kind == 'table':
            markdown += ['| '+ ' | '.join(row)+' |' for row in [body[0], ('---','---'), *body[1:]]]
        else:
            markdown.append(('# ' if kind == 'title' else '## ' if kind == 'heading' else '')+body)
        markdown.append('')
    (ROOT/'report.md').write_text('\n'.join(markdown), encoding='utf-8')
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
    regular = find_font(args.font_regular, ['/System/Library/Fonts/Supplemental/Arial.ttf',
                   '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'])
    bold = find_font(args.font_bold, ['/System/Library/Fonts/Supplemental/Arial Bold.ttf',
                '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'])
    pdfmetrics.registerFont(TTFont('Report', regular))
    pdfmetrics.registerFont(TTFont('ReportBold', bold))
    styles = {
      'title': ParagraphStyle('Title', fontName='ReportBold', fontSize=17, leading=21,
                              textColor=colors.HexColor('#16324f'), spaceAfter=7),
      'subtitle': ParagraphStyle('Subtitle', fontName='Report', fontSize=9, leading=12,
                                 textColor=colors.HexColor('#486174'), spaceAfter=12),
      'heading': ParagraphStyle('Heading', fontName='ReportBold', fontSize=11, leading=14,
                                textColor=colors.HexColor('#16324f'), spaceBefore=8, spaceAfter=5,
                                keepWithNext=True),
      'paragraph': ParagraphStyle('Body', fontName='Report', fontSize=9.4, leading=13.0,
                                  spaceAfter=6, alignment=TA_LEFT),
      'small': ParagraphStyle('Small', fontName='Report', fontSize=7.8, leading=10.3,
                              textColor=colors.HexColor('#33475b'), spaceAfter=6)
    }
    story = []
    for kind, body in blocks:
        if kind == 'pagebreak':
            story.append(PageBreak())
        elif kind == 'table':
            t = Table([[Paragraph(escape(c), styles['small']) for c in row] for row in body],
                      colWidths=[81*mm, 89*mm], hAlign='LEFT')
            t.setStyle(TableStyle([
                ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e3edf4')),
                ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f8fa')]),
                ('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),
                ('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),2),
                ('LINEBELOW',(0,0),(-1,0),0.5,colors.HexColor('#adc2d0'))]))
            story += [t, Spacer(1, 3)]
        else:
            story.append(Paragraph(escape(body).replace('\n','<br/>'), styles[kind]))
    doc = SimpleDocTemplate(str(ROOT/'report.pdf'), pagesize=A4,
             rightMargin=20*mm, leftMargin=20*mm, topMargin=15*mm, bottomMargin=17*mm,
             title='Домашнее задание 4 - вариант 6', author='Данил Белов')
    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont('Report', 8)
        canvas.setFillColor(colors.HexColor('#486174'))
        canvas.drawString(20*mm, 10*mm, 'Биоинформатика | Вариант 6 | GRCh38')
        canvas.drawRightString(A4[0]-20*mm, 10*mm, str(doc.page))
        canvas.restoreState()
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print('Created report.md and report.pdf')


if __name__ == '__main__':
    main()
