#!/bin/bash
# Modified POLCA polishing script for Flyest (MaSuRCA-v4.1.0 derived)
# Fixes include proper read parsing, robust I/O checks, and safer success flags.

MYPATH="$(dirname "$0")"
MYPATH="$(cd "$MYPATH" && pwd)"
UFASTA="$MYPATH/../binaries/ufasta"
FIX_CONSENSUS="$MYPATH/fix_consensus_from_vcf.pl"
set -o pipefail

export NUM_THREADS=4
export MEM=1G
export FIX=1
export BATCH_SIZE=5000000

GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
    log "Aborted"
    kill -9 0
    exit 1
}

log() {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}

error_exit() {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}

usage() {
    echo "Modified polca.sh for Flyest"
    echo "Usage:"
    echo "polca_mod.sh -a <assembly.fasta> -r <reads1.fq reads2.fq> [-t threads] [-m mem] [-n nofix]"
    echo ""
    echo "Requires: bwa, samtools, freebayes, ufasta"
}

# ---------- Parse arguments ----------
if [[ $# -eq 0 ]]; then
    usage; exit 1
fi

READS=()
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -a|--assembly)
            ASM="$2"; shift ;;
        -r|--reads)
            shift
            while [[ $# -gt 0 && $1 != -* ]]; do
                READS+=("$1")
                shift
            done ;;
        -t|--threads)
            NUM_THREADS="$2"; shift ;;
        -m|--memory)
            MEM="$2"; shift ;;
        -n|--nofix)
            FIX=0 ;;
        -b|--batch)
            BATCH_SIZE="$2"; shift ;;
        -v|--verbose)
            set -x ;;
        -h|--help|--usage)
            usage; exit 0 ;;
        *)
            echo "Unknown option $1"; exit 1 ;;
    esac
    shift
done

[[ ! -e $ASM ]] && error_exit "Assembly $ASM not found!"

which bwa >/dev/null || error_exit "bwa not found!"
which freebayes >/dev/null || error_exit "freebayes not found!"
which samtools >/dev/null || error_exit "samtools not found!"
[[ -x "$UFASTA" ]] || error_exit "Bundled ufasta not found or not executable at $UFASTA"
[[ -f "$FIX_CONSENSUS" ]] || error_exit "Missing fix_consensus_from_vcf.pl in $MYPATH"

ASMPATH="$(dirname "$ASM")"
ASMPATH="$(cd "$ASMPATH" && pwd)"
BASM="$(basename "$ASM")"

BWA=$(which bwa)
FREEBAYES=freebayes
SAMTOOLS=samtools

# ---------- Index ----------
if [ ! -e "$ASMPATH/$BASM.index.success" ]; then
    log "Creating BWA index for $ASM"
    rm -f "$ASMPATH/$BASM.map.success"
    $BWA index "$ASM" -p "$ASMPATH/$BASM.bwa" 2>>"$ASMPATH/bwa.err" \
        && touch "$ASMPATH/$BASM.index.success" \
        || error_exit "BWA index failed for $ASM"
fi

# ---------- Mapping ----------
if [ ! -e "$ASMPATH/$BASM.map.success" ]; then
    log "Aligning reads to $ASM"
    rm -f "$ASMPATH/$BASM.sort.success"
    $BWA mem -SP -t $NUM_THREADS "$ASMPATH/$BASM.bwa" "${READS[@]}" \
        > "$ASMPATH/$BASM.unSorted.sam" 2>>"$ASMPATH/bwa.err" \
        || error_exit "BWA mem failed for $BASM"

    [ ! -s "$ASMPATH/$BASM.unSorted.sam" ] && error_exit "Empty SAM file for $BASM"
    touch "$ASMPATH/$BASM.map.success"
fi

# ---------- Sort & index ----------
if [ ! -e "$ASMPATH/$BASM.sort.success" ]; then
    log "Sorting and indexing alignment file"
    rm -f "$ASMPATH/$BASM.vc.success"
    $SAMTOOLS sort -m $MEM -@ $NUM_THREADS -o "$ASMPATH/$BASM.alignSorted.bam" <(samtools view -uhS "$ASMPATH/$BASM.unSorted.sam") 2>>"$ASMPATH/samtools.err" \
        && $SAMTOOLS index "$ASMPATH/$BASM.alignSorted.bam" 2>>"$ASMPATH/samtools.err" \
        && $SAMTOOLS faidx "$ASM" 2>>"$ASMPATH/samtools.err" \
        && touch "$ASMPATH/$BASM.sort.success" \
        || error_exit "Samtools sort/index failed"
    [ ! -s "$ASMPATH/$BASM.alignSorted.bam" ] && error_exit "Empty BAM file for $BASM"
fi

# ---------- Variant calling ----------
if [ ! -e "$ASMPATH/$BASM.vc.success" ]; then
    rm -f "$ASMPATH/$BASM.report.success" "$ASMPATH/$BASM.fix.success"
    log "Calling variants with FreeBayes"
    mkdir -p "$ASMPATH/$BASM.vc"

    "$UFASTA" sizes -H "$ASM" | awk 'BEGIN{bs='$BATCH_SIZE';bi=1;btot=0;}{if($2>bs){for(n=0;n<$2;n+=bs){max=n+bs;if(max>$2) max=$2;print bi" "$1" "n" "max;bi++;}btot=0}else{if(btot+$2>bs){bi++;btot=$2}else{btot+=$2}print bi" "$1" 0 "$2}}' > "$ASMPATH/$BASM.batches"
    BATCHES=$(awk '{print $1}' "$ASMPATH/$BASM.batches" | uniq | wc -l | awk '{print $1}')

    (
        cd "$ASMPATH/$BASM.vc"
        rm -f "$BASM.vc.success" "$BASM.fix.success"
        echo "#!/bin/bash" > commands.sh
        echo "if [ ! -e \$1.vc.success ]; then" >> commands.sh
        echo "awk '{if(\$1=='\$1') print \$2\" \"\$3\" \"\$4}' ../$BASM.batches > batch.\$1" >> commands.sh
        echo "$FREEBAYES -t batch.\$1 -m 0 --min-coverage 3 -R 0 -p 1 -F 0.2 -E 0 -b ../$BASM.alignSorted.bam -v \$1.vcf -f $ASMPATH/$BASM && touch \$1.vc.success" >> commands.sh
        echo "fi" >> commands.sh
        chmod 0755 commands.sh
        seq 1 $BATCHES | xargs -P $NUM_THREADS -I % ./commands.sh %

        for f in $(seq 1 $BATCHES); do
            [ ! -e $f.vc.success ] && error_exit "FreeBayes failed on batch $f"
        done
        touch "$BASM.vc.success"
    )

    if [ -e "$ASMPATH/$BASM.vc/$BASM.vc.success" ]; then
        seq 1 $BATCHES | xargs -I % ls "$ASMPATH/$BASM.vc/%.vcf" | xargs cat | grep '^##' | grep -v commandline | sort -S 10% | uniq > "$ASMPATH/$BASM.vcf.header1"
        seq 1 $BATCHES | xargs -I % ls "$ASMPATH/$BASM.vc/%.vcf" | xargs cat | grep '^#C' | sort -S 10% | uniq > "$ASMPATH/$BASM.vcf.header2"
        seq 1 $BATCHES | xargs -I % ls "$ASMPATH/$BASM.vc/%.vcf" | xargs cat | grep -v '^#' > "$ASMPATH/$BASM.vcf.body"
        cat "$ASMPATH/$BASM.vcf.header1" "$ASMPATH/$BASM.vcf.header2" "$ASMPATH/$BASM.vcf.body" > "$ASMPATH/$BASM.vcf"
        touch "$ASMPATH/$BASM.vc.success"
        rm -rf "$ASMPATH/$BASM.vc"
    else
        error_exit "FreeBayes failed to produce VCF"
    fi
fi

# ---------- Fix consensus ----------
if [ ! -e "$ASMPATH/$BASM.fix.success" ] && [ $FIX -gt 0 ]; then
    mkdir -p "$ASMPATH/$BASM.fix"
    "$UFASTA" sizes -H "$ASM" | sort -S 10% -k2 | awk '{print $1}' > "$ASMPATH/$BASM.names"
    (
        cd "$ASMPATH/$BASM.fix"
        CONTIGS=$(wc -l "$ASMPATH/$BASM.names" | awk '{print $1}')
        BATCH_SIZE=$((CONTIGS / NUM_THREADS + 1))
        [ $BATCH_SIZE -lt 1 ] && BATCH_SIZE=1
        echo "#!/bin/bash" > commands.sh
        echo "if [ ! -e \$1.fix.success ]; then" >> commands.sh
        echo "$FIX_CONSENSUS <($UFASTA extract -f $ASMPATH/$BASM.names $ASMPATH/$BASM) < ../$BASM.vcf > \$1.fixed.tmp 2>\$1.fixed.err && mv \$1.fixed.tmp \$1.fixed && touch \$1.fix.success" >> commands.sh
        echo "fi" >> commands.sh
        chmod 0755 commands.sh
        seq 1 $BATCH_SIZE | xargs -P $NUM_THREADS -I % ./commands.sh %
        for f in $(seq 1 $BATCH_SIZE); do
            [ ! -e $f.fix.success ] && error_exit "Fix failed on batch $f"
        done
        touch "$ASMPATH/$BASM.fix.success"
    )
    [ -e "$ASMPATH/$BASM.fix.success" ] && \
        cat "$ASMPATH/$BASM.fix/"*.fixed | "$UFASTA" format > "$ASMPATH/$BASM.PolcaCorrected.fa" \
        || error_exit "Consensus fix failed"
fi

# ---------- Reporting ----------
log "Creating report"
NUMSUB=$(grep --text -v '^#' "$ASMPATH/$BASM.vcf" | perl -ane '{if(length($F[3])==1 && length($F[4])==1){ print "$F[9]:1\n";}}' | awk -F ':' 'BEGIN{nerr=0}{if($4==0 && $6>1) nerr+=$NF}END{print nerr}')
NUMIND=$(grep --text -v '^#' "$ASMPATH/$BASM.vcf" | perl -ane '{if(length($F[3])>1 || length($F[4])>1){$nerr=abs(length($F[3])-length($F[4]));print "$F[9]:$nerr\n";}}' | awk -F ':' 'BEGIN{nerr=0}{if($4==0 && $6>1) nerr+=$NF}END{print nerr}')
ASMSIZE=$("$UFASTA" n50 -S "$ASM" | awk '{print $2}')
NUMERR=$((NUMSUB + NUMIND))
QUAL=$(echo $NUMERR $ASMSIZE | awk '{print 100-$1/$2*100}')
QV=$(perl -e '{$erate='$NUMERR'/'$ASMSIZE';printf("%.2f\n",-10*log($erate+0.0000000001)/log(10))}')
{
echo "Stats BEFORE polishing:"
echo "Substitution Errors Found: $NUMSUB"
echo "Insertion/Deletion Errors Found: $NUMIND"
echo "Assembly Size: $ASMSIZE"
echo "Consensus Quality Before Polishing: $QUAL"
echo "Consensus QV Before Polishing: $QV"
} > "$ASMPATH/$BASM.report"

# ---------- Cleanup ----------
mkdir -p "$ASMPATH/logs" "$ASMPATH/tmp"
mv "$ASMPATH/$BASM.PolcaCorrected.fa" "$ASMPATH/tmp/${BASM%%.fasta}_polca.fasta" 2>/dev/null
mv "$ASMPATH/$BASM.report" "$ASMPATH/logs/${BASM%%.fasta}_polca.report" 2>/dev/null
rm -rf "$ASMPATH"/*err "$ASMPATH"/*fix "$ASMPATH"/*names "$ASMPATH"/*success "$ASMPATH"/*alignSorted.bam* "$ASMPATH"/*unSorted.sam* "$ASMPATH"/*fasta.vcf "$ASMPATH"/*batches

log "✅ Success! Final report: $ASMPATH/logs/${BASM%%.fasta}_polca.report"
log "✅ Polished assembly: $ASMPATH/tmp/${BASM%%.fasta}_polca.fasta"
